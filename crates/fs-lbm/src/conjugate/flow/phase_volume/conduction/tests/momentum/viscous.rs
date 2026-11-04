//! G1/G3/G4: original conductive heat, physical stress heat and one inventory.

use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::{MaterialMomentumViscosityInput, MaterialViscosityControls};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

fn viscous_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 200,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-8,
        maximum_angular_impulse_defect_kg_m2_s: 1e-8,
        maximum_energy_defect_j: 1e-8,
    }
}

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
fn original_heat_and_viscous_impulse_share_one_pressure_and_material_transport() {
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
    let mut volume = input(&mass, &states, &[0.0; 2]);
    volume.boundary_inflow[0] = Some(caloric.state_at_specific_enthalpy(141.0 / 100.0).unwrap());
    let mut incoming = [None; 6];
    incoming[0] = Some([4.0 / 5.0, 0.0, 0.0]);
    let data = ThermalVolumeMomentumInput {
        volume,
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: incoming,
    };
    let gate = CancelGate::new_clock_free();
    let run = |data| {
        field
        .advance_isotropic_conductive_thermal_volume_momentum_with_viscosity_and_pressure_profile(
            &domain, faces, profile, &conductivity, data, viscosity,
            coupled_controls(1.0 / 8.0), viscous_controls(), &gate,
        )
    };
    let result = run(data).unwrap();
    assert_eq!(molecular.calls.get(), 2);
    assert_eq!(
        result.report.conductivity_source_identity,
        conductivity.identity()
    );
    assert_eq!(result.report.conductive_faces, 1);
    assert_eq!(result.conductive_enthalpy_increment_j, [0.5, -0.5]);
    close(result.report.exchanged_heat_j, 0.5);
    close(result.report.maximum_enthalpy_mixing_weight, 1.0 / 8.0);
    let viscous = &result.evolution.viscosity;
    // Two side walls give stress stiffness 4 per full x face, including
    // the same half weighting at either mouth. rho=2 and dt=1/8 therefore
    // give u_mu=1/(1+4*dt/rho)=4/5, and Q_mu=dt*4*u_mu^2=8/25 per cell.
    for heat in &viscous.cell_molecular_heat_j {
        close(*heat, 8.0 / 25.0);
    }
    close(viscous.molecular_heat_j, 16.0 / 25.0);
    close(viscous.face_impulse_kg_m_s[0], -4.0 / 5.0);
    close(viscous.cell_momentum_change_kg_m_s[0], -4.0 / 5.0);
    close(viscous.backward_euler_loss_j, 2.0 / 25.0);
    close(viscous.numerical_scatter_loss_j, 0.0);
    close(viscous.wall_work_j, 0.0);
    let evolution = &result.evolution.evolution;
    let volume = &evolution.thermal.volume;
    let transport = &volume.transport;
    // Original-state conduction gives +/-1/2 J, before molecular heat.
    // h*=[141/100,291/100]. Each accepted boundary/interior route is 1/5 kg:
    // (2+1/5)h0=2*(141/100)+(1/5)*(141/100),
    // (2+1/5)h1=2*(291/100)+(1/5)*h0.
    for c in 0..2 {
        close(transport.mass_kg[c], 2.0);
        close(
            transport.specific_enthalpy_j_kg[c],
            [141.0 / 100.0, 3051.0 / 1100.0][c],
        );
        close(evolution.momentum.momentum_kg_m_s[c][0], 8.0 / 5.0);
        close(evolution.momentum.velocity_m_s[c][0], 4.0 / 5.0);
        close(volume.volume_source_m3_s[c], 0.0);
        let phase = volume.phase_states[c].unwrap();
        assert_eq!(phase.phase_curve_identity(), caloric.identity());
        assert_eq!(phase.material_card_identity(), material.content_hash());
        close(transport.mass_kg[c] / phase.bulk_density_kg_m3(), 1.0);
    }
    for face in 0..3 {
        close(transport.projection.field.flux_x(face, 0, 0), 4.0 / 5.0);
    }
    close(transport.report.incoming_mass_kg, 1.0 / 5.0);
    close(transport.report.outgoing_mass_kg, 1.0 / 5.0);
    assert_eq!(transport.report.mass_routes, 3);
    close(volume.report.requested_enthalpy_increment_j, 16.0 / 25.0);
    close(
        volume.report.material_enthalpy_change_j,
        16.0 / 25.0 + (141.0 / 100.0 - 3051.0 / 1100.0) / 5.0,
    );
    close(evolution.force_projection.pressure_impulse_kg_m_s[0], 0.0);
    close(evolution.feedback.cell_momentum_change_kg_m_s[0], 0.0);
    assert_eq!(
        evolution.momentum.report.represented_impulse_kg_m_s,
        [0.0; 3]
    );
    close(evolution.momentum.report.numerical_transport_loss_j, 0.0);
    // Physical molecular heat is deposited once; the BE loss remains separate.
    close(
        evolution.momentum.report.kinetic_energy_after_j - 2.0
            + viscous.molecular_heat_j
            + viscous.backward_euler_loss_j,
        0.0,
    );
    close(result.report.absolute_enthalpy_defect_j, 0.0);

    let missing_velocity = ThermalVolumeMomentumInput {
        boundary_velocity_m_s: [None; 6],
        ..data
    };
    assert!(
        matches!(run(missing_velocity), Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("incoming mass route lacks an authored velocity"))
    );
    assert_eq!(molecular.calls.get(), 4);
    assert_eq!(run(data).unwrap(), result);
    assert_eq!(molecular.calls.get(), 6);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    assert_eq!(
        states,
        [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()))
    );
}

#[test]
fn conductive_and_molecular_addition_defects_share_one_total_enthalpy_budget() {
    let span = 2.0_f64.powi(60);
    let authored = 2.0_f64.powi(55);
    let material = card("T", &[(300.0, 4.0), (304.0, 4.0)]);
    let caloric = curve(
        &material,
        &[(0.0, 300.0, 1.0, 2.0), (span, 304.0, 1.0, 2.0)],
    );
    let conductivity = source(&material, &caloric, [0.0, span]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let states = [0.0, span / 4.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let momentum = [[2.0, 0.0, 0.0]; 2];
    let force = [[0.0; 3]; 2];
    let authored_heat = [authored; 2];
    let data = ThermalVolumeMomentumInput {
        volume: input(&mass, &states, &authored_heat),
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: [None; 6],
    };
    let molecular = FluidTransportState::try_new(2.0, 1.0, 4.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&molecular); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Symmetry; 6],
    };
    let mut viscous = viscous_controls();
    viscous.maximum_energy_defect_j = 2.0;
    let mut loose = coupled_controls(1.0 / 8.0);
    loose
        .evolution
        .volume
        .thermal
        .volume
        .transport
        .maximum_enthalpy_defect_j = 2.0;
    let gate = CancelGate::new_clock_free();
    let run = |controls| {
        field.advance_isotropic_conductive_thermal_volume_momentum_with_viscosity(
            &domain,
            &conductivity,
            data,
            viscosity,
            controls,
            viscous,
            &gate,
        )
    };
    let result = run(loose).unwrap();
    // T=[300,301], G=4 and dt=1/8 give paired heat +/-1/2. Both
    // disappear beside authored 2^55: their absolute addition defect is 1.
    // The sole interior x face has M=2, stress stiffness 4 and u_mu=4/5.
    // Its two normal strains each dissipate dt*2*u_mu^2=4/25 J; those
    // also disappear beside 2^55. Subsequent binary-power calorics are exact.
    assert_eq!(result.conductive_enthalpy_increment_j, [0.5, -0.5]);
    close(result.report.increment_assembly_defect_j, 1.0);
    close(result.evolution.viscosity.molecular_heat_j, 8.0 / 25.0);
    close(result.evolution.absolute_heat_assembly_defect_j, 8.0 / 25.0);
    close(
        result
            .evolution
            .evolution
            .thermal
            .volume
            .report
            .absolute_enthalpy_defect_j,
        0.0,
    );
    close(result.report.absolute_enthalpy_defect_j, 1.0 + 8.0 / 25.0);
    let mut tight = loose;
    tight
        .evolution
        .volume
        .thermal
        .volume
        .transport
        .maximum_enthalpy_defect_j = 1.1;
    let conduction_only = field
        .advance_isotropic_conductive_thermal_volume_momentum(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            data,
            tight,
            &gate,
        )
        .unwrap();
    close(conduction_only.report.absolute_enthalpy_defect_j, 1.0);
    let viscosity_only = field
        .advance_piecewise_thermal_volume_momentum_with_viscosity(
            &domain,
            &caloric,
            data,
            viscosity,
            tight.evolution,
            viscous,
            &gate,
        )
        .unwrap();
    close(viscosity_only.absolute_enthalpy_defect_j, 8.0 / 25.0);
    assert!(
        matches!(run(tight), Err(ChtError::InvalidInput { field: "flow.conductive_thermal_volume", reason })
        if reason.contains("conduction and native volume transport exceed"))
    );
    assert_eq!(run(loose).unwrap(), result);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
}

#[test]
fn incoming_material_cannot_hide_molecular_heating_outside_conductivity_support() {
    let material = card("T", &[(300.0, 1.0), (302.0, 1.0)]);
    // v(h)=1/2+h/4; molecular heat changes physical density and volume.
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (2.0, 302.0, 1.0, 1.0)]);
    let narrow = source(&material, &caloric, [0.0, 0.23]);
    let broad = source(&material, &caloric, [0.0, 0.3]);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0];
    let states = [Some(caloric.state_at_specific_enthalpy(0.0).unwrap())];
    let momentum = [[2.0, 0.0, 0.0]];
    let force = [[0.0; 3]];
    let molecular = FluidTransportState::try_new(2.0, 1.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 1] = [Some(&molecular)];
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
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let profile = StaggeredBoundaryPressure {
        x_pa: &[Some(0.0), Some(0.0)],
        y_pa: &[None; 2],
        z_pa: &[None; 2],
    };
    let mut volume = input(&mass, &states, &[0.0]);
    volume.boundary_inflow[0] = states[0];
    let mut incoming = [None; 6];
    incoming[0] = Some([0.5, 0.0, 0.0]);
    let data = ThermalVolumeMomentumInput {
        volume,
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &force,
        boundary_velocity_m_s: incoming,
    };
    let gate = CancelGate::new_clock_free();
    let run = |conductivity| {
        field
        .advance_isotropic_conductive_thermal_volume_momentum_with_viscosity_and_pressure_profile(
            &domain, faces, profile, conductivity, data, viscosity,
            coupled_controls(0.5), viscous_controls(), &gate,
        )
    };
    let accepted = run(&broad).unwrap();
    assert_eq!(accepted.report.conductive_faces, 0);
    close(accepted.report.exchanged_heat_j, 0.0);
    close(accepted.evolution.viscosity.molecular_heat_j, 0.5);
    // Viscous speed is 1/2; the deposited .5 J reaches h*=.25 before
    // transport. Expansion supplies S=.25. With Q_left=3/8 and Q_right=5/8,
    // inlet mass is 3/8 and M'=(2+3/8)/(1+5/16)=38/21. Incoming h=0
    // cools the endpoint to h'=(.5)/(2+3/8)=4/19, inside the narrow source.
    let transport = &accepted.evolution.evolution.thermal.volume.transport;
    close(transport.mass_kg[0], 38.0 / 21.0);
    close(transport.specific_enthalpy_j_kg[0], 4.0 / 19.0);
    assert!(transport.specific_enthalpy_j_kg[0] < narrow.enthalpy_bounds_j_kg()[1]);
    assert!(
        matches!(run(&narrow), Err(ChtError::InvalidInput { reason, .. })
        if reason.contains("molecular heating leaves the original conductivity source"))
    );
    assert_eq!(run(&broad).unwrap(), accepted);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]]);
    assert_eq!(
        states,
        [Some(caloric.state_at_specific_enthalpy(0.0).unwrap())]
    );
}
