//! G1/G3: finite annular thermal walls, independent radial equations and work.
use fs_duct::DuctError;
use fs_duct::fluid::FluidDuctError;
use fs_duct::fluid::annular::thermal_wall::AnnularThermalWall as Wall;
use fs_duct::fluid::annular::{AnnularAdmission, AnnularSection};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::{PI, TAU};

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn admission() -> AnnularAdmission {
    AnnularAdmission {
        maximum_transverse_phase: 0.8,
        maximum_cell_diffusion_phase: 0.9,
        maximum_relative_residual: 1e-7,
        maximum_payload_bytes: 1_000_000,
    }
}
fn finite(g: f64, c: f64) -> Wall {
    Wall::Finite {
        conductance_w_m2_k: g,
        areal_heat_capacity_j_m2_k: c,
    }
}
fn relative(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}
fn near(actual: C64, expected: C64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * expected.abs().max(1e-300),
        "{actual:?} != {expected:?}"
    );
}
fn material(count: usize) -> RelaxationState {
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
fn isothermal_limit_recovers_the_original_annular_wave_and_empty_reservoir_parity() {
    let material = material(0);
    let section = AnnularSection::new(0.002, 0.003, 256).unwrap();
    for omega in [20.0, 500.0, 5000.0] {
        let old = section
            .evaluate(material.gas(), omega, admission())
            .unwrap();
        let fields = section
            .prepare_thermal_wall_fields(material.gas(), omega, [Wall::Isothermal; 2], admission())
            .unwrap();
        near(fields.wave().wavenumber, old.wave.wavenumber, 3e-12);
        near(
            fields.wave().characteristic_impedance,
            old.wave.characteristic_impedance,
            3e-12,
        );
        near(
            fields.mean_temperature_response(),
            old.thermal.unwrap().mean,
            3e-12,
        );
        let molecular = section
            .prepare_molecular_thermal_wall_fields(
                &material,
                omega,
                [Wall::Isothermal; 2],
                admission(),
            )
            .unwrap();
        assert_eq!(fields.wave().wavenumber, molecular.wave().wavenumber);
        assert!(std::ptr::eq(molecular.source(), &material));
        for side in 0..2 {
            let wall = fields.wall_observation(side, C64::ONE).unwrap();
            assert_eq!(wall.temperature_k, C64::ZERO);
            assert_eq!(wall.thermal_storage_j_m, 0.0);
            assert_eq!(wall.bath_dissipation_w_m, 0.0);
            assert_eq!(
                wall.fluid_to_wall_heat_flux_w_m2,
                wall.wall_to_bath_heat_flux_w_m2
            );
            assert!(wall.fluid_to_wall_heat_flux_w_m2.abs() > 0.0);
        }
    }
}

#[test]
fn insulated_walls_preserve_uniform_temperature_and_shear_without_thermal_loss() {
    let gas = gas();
    let section = AnnularSection::new(0.002, 0.003, 96).unwrap();
    let fields = section
        .prepare_thermal_wall_fields(&gas, 300.0, [Wall::Adiabatic; 2], admission())
        .unwrap();
    let p = C64::new(3.0, -1.0);
    let expected = p / C64::from_re(gas.density * gas.specific_heat_cp);
    for cell in [0, 17, 60, 95] {
        let value = fields
            .field_in_cell(cell, 0.31, p, C64::new(1e-8, -2e-8))
            .unwrap();
        near(value.gas_temperature_k, expected, 4e-15);
        assert_eq!(value.radial_temperature_gradient_k_m, C64::ZERO);
        assert_eq!(value.conduction_dissipation_w_m3, 0.0);
    }
    let work = fields
        .dissipation_per_axial_length(p, C64::from_re(1e-8))
        .unwrap();
    assert!(work.viscous_w_m > 0.0);
    assert_eq!(work.conduction_w_m, 0.0);
    assert_eq!(work.molecular_w_m, 0.0);
    assert_eq!(work.walls_w_m, [0.0; 2]);
    for side in 0..2 {
        let wall = fields.wall_observation(side, p).unwrap();
        assert_eq!(wall.fluid_to_wall_heat_flux_w_m2, C64::ZERO);
        assert_eq!(wall.wall_to_bath_heat_flux_w_m2, C64::ZERO);
        assert_eq!(wall.thermal_storage_j_m, 0.0);
    }
}

// Independent dimensional ODE shooting, with fixed-step RK4. Production uses
// weighted P1 weak equations and a tridiagonal solve. Third state integrates
// the physical annular area, so no production mass/stiffness entry is reused.
fn admittance(wall: Wall, omega: f64) -> C64 {
    match wall {
        Wall::Finite {
            conductance_w_m2_k: g,
            areal_heat_capacity_j_m2_k: c,
        } => C64::new(g, omega * c),
        _ => C64::ZERO,
    }
}
fn ode_solution(gas: &GasState, cp: C64, omega: f64, a: f64, b: f64, walls: [Wall; 2]) -> [C64; 3] {
    let z = C64::new(cp.im, cp.re).scale(omega * gas.density / gas.thermal_conductivity);
    let yi = admittance(walls[0], omega).scale(1.0 / gas.thermal_conductivity);
    let yo = admittance(walls[1], omega).scale(1.0 / gas.thermal_conductivity);
    let basis = if walls[0] == Wall::Isothermal {
        [C64::ZERO, C64::ONE, C64::ZERO]
    } else {
        [C64::ONE, yi, C64::ZERO]
    };
    let integrate = |mut v: [C64; 3], forcing: bool| {
        let steps = 4096;
        let h = (b - a) / steps as f64;
        let rhs = |r: f64, u: [C64; 3]| {
            [
                u[1],
                z * (u[0] - if forcing { C64::ONE } else { C64::ZERO }) - u[1].scale(1.0 / r),
                u[0].scale(TAU * r),
            ]
        };
        let advance = |u: [C64; 3], slope: [C64; 3], dt: f64| {
            core::array::from_fn(|i| u[i] + slope[i].scale(dt))
        };
        for step in 0..steps {
            let r = a + h * step as f64;
            let k1 = rhs(r, v);
            let k2 = rhs(r + 0.5 * h, advance(v, k1, 0.5 * h));
            let k3 = rhs(r + 0.5 * h, advance(v, k2, 0.5 * h));
            let k4 = rhs(r + h, advance(v, k3, h));
            v = core::array::from_fn(|i| {
                v[i] + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]).scale(h / 6.0)
            });
        }
        v
    };
    let drive = integrate([C64::ZERO; 3], true);
    let free = integrate(basis, false);
    let residual = |v: [C64; 3]| {
        if walls[1] == Wall::Isothermal {
            v[0]
        } else {
            v[1] + yo * v[0]
        }
    };
    let amplitude = -residual(drive) / residual(free);
    [
        (drive[2] + amplitude * free[2]).scale(1.0 / (PI * (b * b - a * a))),
        amplitude * basis[0],
        drive[0] + amplitude * free[0],
    ]
}

#[test]
fn finite_and_mixed_walls_converge_to_independent_dimensional_radial_equations() {
    let material = material(2);
    let gas = material.gas();
    let omega = 400.0;
    let cp = material.caloric_response(omega).unwrap().cp_j_kg_k;
    for walls in [
        [finite(100.0, 0.05), finite(12.0, 0.2)],
        [Wall::Isothermal, finite(30.0, 0.08)],
        [Wall::Adiabatic, Wall::Isothermal],
    ] {
        let expected = ode_solution(gas, cp, omega, 0.002, 0.003, walls);
        let mut previous = f64::INFINITY;
        for elements in [16, 32, 64] {
            let section = AnnularSection::new(0.002, 0.003, elements).unwrap();
            let fields = section
                .prepare_molecular_thermal_wall_fields(&material, omega, walls, admission())
                .unwrap();
            let error = (fields.mean_temperature_response() - expected[0]).abs();
            assert!(
                error < 0.3 * previous,
                "{walls:?}: {elements}: {error:e} >= {previous:e}"
            );
            previous = error;
            if elements == 64 {
                near(fields.mean_temperature_response(), expected[0], 3e-4);
                for side in 0..2 {
                    let expected_t = expected[side + 1].conj() / cp / C64::from_re(gas.density);
                    let actual = fields
                        .wall_observation(side, C64::ONE)
                        .unwrap()
                        .temperature_k;
                    if walls[side] == Wall::Isothermal {
                        assert_eq!(actual, C64::ZERO);
                    } else {
                        near(actual, expected_t, 5e-4);
                    }
                }
            }
        }
    }
}

#[test]
fn wall_storage_heat_flux_and_each_positive_loss_close_the_original_native_power_balance() {
    let material = material(2);
    let section = AnnularSection::new(0.002, 0.003, 192).unwrap();
    let walls = [finite(80.0, 0.1), finite(15.0, 0.3)];
    let (p, q) = (C64::new(3.0, -1.0), C64::new(2e-8, 3e-8));
    for omega in [20.0, 400.0, 5000.0] {
        let fields = section
            .prepare_molecular_thermal_wall_fields(&material, omega, walls, admission())
            .unwrap();
        let work = fields.dissipation_per_axial_length(p, q).unwrap();
        let wave = fields.wave();
        let reference = 0.5 * (wave.wavenumber * wave.characteristic_impedance).im * q.norm_sq()
            + 0.5 * (wave.wavenumber / wave.characteristic_impedance).im * p.norm_sq();
        relative(work.total_w_m, reference, 2e-10);
        assert!(work.viscous_w_m > 0.0 && work.conduction_w_m > 0.0 && work.molecular_w_m > 0.0);
        for side in 0..2 {
            let observation = fields.wall_observation(side, p).unwrap();
            let Wall::Finite {
                conductance_w_m2_k: g,
                areal_heat_capacity_j_m2_k: c,
            } = walls[side]
            else {
                unreachable!()
            };
            near(
                observation.fluid_to_wall_heat_flux_w_m2,
                observation.temperature_k * C64::new(g, -omega * c),
                3e-15,
            );
            near(
                observation.wall_to_bath_heat_flux_w_m2,
                observation.temperature_k.scale(g),
                3e-15,
            );
            relative(
                observation.bath_dissipation_w_m,
                PI * observation.radius_m * g / material.gas().temperature
                    * observation.temperature_k.norm_sq(),
                3e-15,
            );
            relative(
                observation.thermal_storage_j_m,
                PI * observation.radius_m * c / (2.0 * material.gas().temperature)
                    * observation.temperature_k.norm_sq(),
                3e-15,
            );
            assert_eq!(observation.bath_dissipation_w_m, work.walls_w_m[side]);
            assert!(observation.thermal_storage_j_m > 0.0);
        }
    }
}

#[test]
fn local_physical_volume_quadrature_recovers_all_three_fluid_work_channels() {
    let material = material(2);
    let section = AnnularSection::new(0.002, 0.003, 128).unwrap();
    let fields = section
        .prepare_molecular_thermal_wall_fields(
            &material,
            600.0,
            [finite(30.0, 0.15), Wall::Adiabatic],
            admission(),
        )
        .unwrap();
    let (p, q) = (C64::new(2.0, -0.3), C64::new(3e-8, -2e-8));
    let mut integral = [0.0; 3];
    for cell in 0..section.elements() {
        let dr = section.radius_m(cell + 1).unwrap() - section.radius_m(cell).unwrap();
        for fraction in [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()] {
            let field = fields.field_in_cell(cell, fraction, p, q).unwrap();
            let weight = PI * field.radius_m * dr;
            for (sum, value) in integral.iter_mut().zip([
                field.viscous_dissipation_w_m3,
                field.conduction_dissipation_w_m3,
                field.molecular_dissipation_w_m3,
            ]) {
                *sum += weight * value;
            }
        }
    }
    let work = fields.dissipation_per_axial_length(p, q).unwrap();
    for (integral, exact) in
        integral
            .into_iter()
            .zip([work.viscous_w_m, work.conduction_w_m, work.molecular_w_m])
    {
        relative(integral, exact, 3e-13);
    }
}

#[test]
fn isolated_finite_wall_capacity_has_the_correct_quasistatic_shared_temperature() {
    let gas = gas();
    let section = AnnularSection::new(0.002, 0.0025, 32).unwrap();
    let walls = [finite(0.0, 0.4), finite(0.0, 0.7)];
    let fluid_capacity = section.area_m2() * gas.density * gas.specific_heat_cp;
    let wall_capacity = TAU * (0.002 * 0.4 + 0.0025 * 0.7);
    let expected = C64::from_re(fluid_capacity / (fluid_capacity + wall_capacity));
    let mut previous = f64::INFINITY;
    for omega in [10.0, 1.0, 0.1] {
        let fields = section
            .prepare_thermal_wall_fields(&gas, omega, walls, admission())
            .unwrap();
        let error = (fields.mean_temperature_response() - expected).abs();
        assert!(
            error < 0.12 * previous,
            "{omega}: {error:e} >= {previous:e}"
        );
        previous = error;
        let work = fields
            .dissipation_per_axial_length(C64::ONE, C64::ZERO)
            .unwrap();
        assert_eq!(work.walls_w_m, [0.0; 2]);
        assert!(
            fields
                .wall_observation(0, C64::ONE)
                .unwrap()
                .thermal_storage_j_m
                > 0.0
        );
        if omega == 0.1 {
            near(fields.mean_temperature_response(), expected, 0.003);
        }
    }
}

#[test]
fn finite_conductance_reaches_insulating_and_isothermal_endpoints_continuously() {
    let gas = gas();
    let section = AnnularSection::new(0.002, 0.003, 192).unwrap();
    let omega = 300.0;
    let insulated = section
        .prepare_thermal_wall_fields(&gas, omega, [Wall::Adiabatic; 2], admission())
        .unwrap();
    let cold = section
        .prepare_thermal_wall_fields(&gas, omega, [Wall::Isothermal; 2], admission())
        .unwrap();
    let almost_insulated = section
        .prepare_thermal_wall_fields(&gas, omega, [finite(1e-4, 0.0); 2], admission())
        .unwrap();
    let almost_cold = section
        .prepare_thermal_wall_fields(&gas, omega, [finite(1e7, 0.0); 2], admission())
        .unwrap();
    near(
        almost_insulated.mean_temperature_response(),
        insulated.mean_temperature_response(),
        2e-6,
    );
    near(
        almost_cold.mean_temperature_response(),
        cold.mean_temperature_response(),
        2e-5,
    );
}

#[test]
fn thermal_boundary_resource_source_and_cancellation_refusals_are_explicit() {
    let gas = gas();
    let section = AnnularSection::new(0.002, 0.003, 128).unwrap();
    for wall in [
        finite(-1.0, 0.0),
        finite(0.0, -1.0),
        finite(f64::NAN, 1.0),
        finite(1.0, f64::INFINITY),
    ] {
        assert!(
            section
                .prepare_thermal_wall_fields(&gas, 200.0, [wall; 2], admission())
                .is_err()
        );
    }
    let mut small = admission();
    small.maximum_payload_bytes = 192 * 129 + 4096;
    assert!(
        section
            .prepare_thermal_wall_fields(&gas, 200.0, [Wall::Adiabatic; 2], small)
            .is_err()
    );
    assert!(matches!(
        section.prepare_thermal_wall_fields_controlled(
            &gas,
            200.0,
            [Wall::Adiabatic; 2],
            admission(),
            || true
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let mut polls = 0;
    assert!(matches!(
        section.prepare_thermal_wall_fields_controlled(
            &gas,
            200.0,
            [finite(30.0, 0.2); 2],
            admission(),
            || {
                polls += 1;
                polls == 8
            }
        ),
        Err(FluidDuctError::Duct(DuctError::Cancelled { .. }))
    ));
    let fields = section
        .prepare_thermal_wall_fields(&gas, 200.0, [finite(30.0, 0.2); 2], admission())
        .unwrap();
    assert!(std::ptr::eq(fields.source(), &gas));
    assert!(std::ptr::eq(fields.section(), &section));
    assert!(fields.field_in_cell(128, 0.0, C64::ONE, C64::ZERO).is_err());
    assert!(
        fields
            .field_in_cell(0, f64::NAN, C64::ONE, C64::ZERO)
            .is_err()
    );
    assert!(fields.wall_observation(2, C64::ONE).is_err());
    assert!(
        fields
            .dissipation_per_axial_length(C64::new(f64::NAN, 0.0), C64::ZERO)
            .is_err()
    );
}
