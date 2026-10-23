//! G1 cross-owner regression: source -> FEEC field -> existing current/heat/flow.
//! All coefficients are synthetic. This is one-way, imposed-field coupling in
//! the low-induction regime, not magnetic-induction or total-energy validation.
mod common;
use common::*;
use fs_em::{
    planar::PlanarMagnetostatics,
    transfer::{CellSample, PlanarFrame},
};
use fs_exec::CancelGate;
use fs_lbm::{
    conjugate::{ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain},
    material::electrothermal::{ElectricHeatInput, ElectricHeatLimits},
    material::{MaterialSlab3, MaterialThermalSlab3, plan_scaling_from_material},
    thermal::sources::HeatSourceLimits,
};
use fs_material::fluid_transport::FluidTransportState;

const DX: f64 = 0.01;
fn fluid() -> MaterialSlab3 {
    let material = FluidTransportState::try_new(1000.0, 1000.0, 100.0, 1.0).unwrap();
    let mut slab = plan_scaling_from_material(&material, DX, 0.01, 0.1)
        .unwrap()
        .thermal_slab3(MaterialThermalSlab3 {
            nx: 2,
            ny: 1,
            nz: 6,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 0.0,
            expansion_per_k: 0.0,
            max_cells: 16,
        })
        .unwrap();
    for cell in 2..14 {
        slab.lattice.grid.f[cell] = fs_lbm::equilibrium3(1.0, [0.02, 0.0, 0.0]);
    }
    slab
}
fn problem(left_mu: f64) -> PlanarMagnetostatics {
    PlanarMagnetostatics::new(
        layered(DX, DX, 8.0 * DX),
        vec![
            material(left_mu),
            material(left_mu),
            material(1e-6),
            material(1e-6),
        ],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn permeability_data_changes_real_current_joule_heat_and_fluid_evolution() {
    let reference = problem(1e-6);
    let changed = problem(2e-6); // This is the only changed physical input.
    let boundary = ends(0.04);
    let original = reference
        .solve(&[0.0; 4], &boundary, config(), &mut || false)
        .unwrap();
    let swapped = changed
        .solve(&[0.0; 4], &boundary, config(), &mut || false)
        .unwrap();
    let domain = VoxelDomain::new(2, 1, 8, DX).unwrap(); // Same dimensions, including slab wall planes.
    let frame = PlanarFrame::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap();
    let samples: Vec<_> = (0..16)
        .map(|cell| {
            let [x, y, z] = domain.coords(cell);
            CellSample {
                triangle: 2 * x,
                position_m: domain.center(x, y, z),
            }
        })
        .collect();
    let a = original
        .sample_cells(&frame, &samples, 1e-12, 16, &mut || false)
        .unwrap();
    let b = swapped
        .sample_cells(&frame, &samples, 1e-12, 16, &mut || false)
        .unwrap();
    for cell in 0..16 {
        near(a[cell][1], -2.0);
        near(
            b[cell][1],
            if cell % 2 == 0 {
                -8.0 / 3.0
            } else {
                -4.0 / 3.0
            },
        );
        near(a[cell][0], 0.0);
        near(b[cell][0], 0.0);
    }
    // Tangential B may jump at the source-material interface. Here B is
    // parallel to that interface, so these sampled slab fields also have zero
    // discrete normal flux divergence; no general conservative remap is claimed.
    let electric = ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 1],
                hi: [2, 1, 2],
                drive: ElectrodeDrive::Voltage(2.0),
            },
            Electrode {
                lo: [0, 0, 6],
                hi: [2, 1, 7],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    };
    let input = ElectricHeatInput {
        setup: &electric,
        contacts: &[],
        limits: ElectricHeatLimits {
            thermal: HeatSourceLimits {
                max_cells: 16,
                temperature_range_k: [290.0, 320.0],
                balance_tolerance_j: 1e-9,
            },
            coupling_tolerance_j: 1e-10,
        },
    };
    let gate = CancelGate::new();
    let mut x = fluid();
    let mut y = fluid();
    for _ in 0..2 {
        let rx = x
            .step_magnetothermal(input, &a, &[[0.0; 3]; 16], &gate, |_, _, _| Ok(Some(1.0)))
            .unwrap();
        let ry = y
            .step_magnetothermal(input, &b, &[[0.0; 3]; 16], &gate, |_, _, _| Ok(Some(1.0)))
            .unwrap();
        assert_ne!(rx.motion.electric.joule_w, ry.motion.electric.joule_w);
        assert_ne!(
            rx.motion.lorentz_force_density_n_m3,
            ry.motion.lorentz_force_density_n_m3
        );
        near(rx.joule_energy_j, rx.heat.source_j);
        near(ry.joule_energy_j, ry.heat.source_j);
        assert!(rx.motion.power_defect_w.abs() < 1e-8);
        assert!(ry.motion.power_defect_w.abs() < 1e-8);
    }
    assert_ne!(x.lattice.grid.f, y.lattice.grid.f); // Actual momentum populations.
    assert_ne!(x.lattice.gpop, y.lattice.gpop); // Actual thermal populations.
    assert_ne!(x.lattice.grid.moments3(8).1, y.lattice.grid.moments3(8).1);
    assert_ne!(
        x.lattice.temperature(0, 0, 3),
        y.lattice.temperature(0, 0, 3)
    );
}
