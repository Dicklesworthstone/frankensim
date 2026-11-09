//! Public material/field/port integration; G1/G3 acoustic work and wall limits.
use super::*;
use fs_material::gas::GasSpec;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> RectangularAdmission {
    RectangularAdmission {
        maximum_relative_tail: 2e-5,
        maximum_transverse_phase: 0.8,
    }
}
fn limits() -> RectangularThermalLimits {
    RectangularThermalLimits {
        cells: [24, 32],
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-7,
        maximum_payload_bytes: 512 * 1024 * 1024,
        sparse: DirectLimits {
            max_symbolic_entries: 8_000_000,
            max_analysis_work: 500_000_000,
            max_numeric_work: 100_000_000_000,
            max_front_dim: 4096,
        },
    }
}
fn finite_wall(g: f64, c: f64) -> RectangularThermalWall {
    RectangularThermalWall::Finite {
        conductance_w_m2_k: g,
        areal_heat_capacity_j_m2_k: c,
    }
}
fn molecular(count: usize) -> RelaxationState {
    let species = SpeciesId::new("synthetic").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * i as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (0.002 + 0.01 * i as f64),
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas(), vec![(species, 1.0)]).unwrap())
        .unwrap()
}
#[test]
fn original_material_four_wall_fields_close_separate_physical_port_losses() {
    let material = molecular(2);
    let section = RectangularSection::new(0.002, 0.003, 256).unwrap();
    let walls = [
        finite_wall(25.0, 0.08),
        finite_wall(80.0, 0.02),
        RectangularThermalWall::Adiabatic,
        RectangularThermalWall::Isothermal,
    ];
    let fields = section
        .prepare_molecular_thermal_wall_fields(&material, 400.0, walls, admission(), limits())
        .unwrap();
    assert!(core::ptr::eq(fields.source(), &material));
    let p = C64::new(1.2, -0.4);
    let q = C64::new(2e-8, 1e-8);
    let heat = fields.dissipation_per_axial_length(p, q).unwrap();
    let wave = fields.wave();
    let series = C64::new(0.0, 1.0) * (wave.wavenumber * wave.characteristic_impedance).conj();
    let shunt = C64::new(0.0, 1.0) * (wave.wavenumber / wave.characteristic_impedance).conj();
    let expected = 0.5 * (series.re * q.norm_sq() + shunt.re * p.norm_sq());
    assert!((heat.total_w_m - expected).abs() < 1e-6 * expected.abs());
    assert!(heat.viscous_w_m > 0.0 && heat.conduction_w_m > 0.0 && heat.molecular_w_m > 0.0);
    assert!(heat.walls_w_m[0] > 0.0 && heat.walls_w_m[1] > 0.0);
    assert_eq!(heat.walls_w_m[2], 0.0);
    assert_eq!(heat.walls_w_m[3], 0.0);
    let zero = fields
        .field_at([0.001, 0.0015], C64::ZERO, C64::ZERO)
        .unwrap();
    assert_eq!(zero.gas_temperature_k, C64::ZERO);
    assert_eq!(zero.axial_velocity_m_s, C64::ZERO);
    for side in 0..4 {
        let w = fields.wall_observation(side, 0.5, p).unwrap();
        if side < 2 {
            assert!(
                (w.fluid_to_wall_heat_flux_w_m2
                    - w.temperature_k * walls[side].coefficients(400.0).unwrap().conj())
                .abs()
                    < 1e-14
            );
            assert!(w.thermal_storage_j_m2 > 0.0);
        }
        if side == 2 {
            assert_eq!(w.fluid_to_wall_heat_flux_w_m2, C64::ZERO);
        }
        if side == 3 {
            assert!(w.temperature_k.abs() < 1e-18);
        }
    }
    assert!(
        fields
            .two_port(0.1, wave.characteristic_impedance.re)
            .is_ok()
    );
}
#[test]
fn empty_reservoirs_and_axis_authorship_preserve_original_gas_fields() {
    let gas = gas();
    let material = molecular(0);
    let a = RectangularSection::new(0.002, 0.003, 256).unwrap();
    let b = RectangularSection::new(0.003, 0.002, 256).unwrap();
    let walls = [
        finite_wall(25.0, 0.08),
        finite_wall(80.0, 0.02),
        RectangularThermalWall::Adiabatic,
        RectangularThermalWall::Isothermal,
    ];
    let x = a
        .prepare_thermal_wall_fields(&gas, 400.0, walls, admission(), limits())
        .unwrap();
    let y = b
        .prepare_molecular_thermal_wall_fields(&material, 400.0, walls, admission(), limits())
        .unwrap();
    assert_eq!(x.mean_temperature_response(), y.mean_temperature_response());
    assert_eq!(x.wave().wavenumber, y.wave().wavenumber);
    assert_eq!(
        x.field_at([0.0007, 0.0023], C64::ONE, C64::new(1e-8, 2e-8))
            .unwrap(),
        y.field_at([0.0007, 0.0023], C64::ONE, C64::new(1e-8, 2e-8))
            .unwrap()
    );
    let insulated = a
        .prepare_thermal_wall_fields(
            &gas,
            400.0,
            [RectangularThermalWall::Adiabatic; 4],
            admission(),
            limits(),
        )
        .unwrap();
    let work = insulated
        .dissipation_per_axial_length(C64::ONE, C64::ZERO)
        .unwrap();
    assert_eq!(work.conduction_w_m, 0.0);
    assert_eq!(work.molecular_w_m, 0.0);
    assert_eq!(work.total_w_m, 0.0);
}
#[test]
fn original_source_and_mesh_admission_refuse_without_prepared_fields() {
    let mut gas = gas();
    let section = RectangularSection::new(0.002, 0.003, 256).unwrap();
    gas.density *= 2.0;
    assert!(
        section
            .prepare_thermal_wall_fields(
                &gas,
                400.0,
                [RectangularThermalWall::Adiabatic; 4],
                admission(),
                limits()
            )
            .is_err()
    );
    let gas = super::tests::gas();
    assert!(
        section
            .prepare_thermal_wall_fields(
                &gas,
                400.0,
                [finite_wall(-1.0, 0.0); 4],
                admission(),
                limits()
            )
            .is_err()
    );
    assert!(
        section
            .prepare_thermal_wall_fields_controlled(
                &gas,
                400.0,
                [RectangularThermalWall::Isothermal; 4],
                admission(),
                limits(),
                || true
            )
            .is_err()
    );
    let mut short = limits();
    short.maximum_payload_bytes = 1;
    assert!(
        section
            .prepare_thermal_wall_fields(
                &gas,
                400.0,
                [RectangularThermalWall::Isothermal; 4],
                admission(),
                short
            )
            .is_err()
    );
}
