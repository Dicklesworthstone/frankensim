use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{
    NsfBoundaryDrive, NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfCharacteristicPatch,
    NsfSgsConfig, NsfSpec, WallThermal,
};
use crate::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_material::gas::{GasSpec, GasState};

fn grid(axis: usize, inlet_high: bool) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let reference = GasState::try_new(&gas, 300.0, 100_000.0).unwrap();
    let unit = (3.0 * reference.specific_gas_constant * 300.0).sqrt();
    let mut dims = [3; 3];
    dims[axis] = 8;
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let mut velocity = [0.0; 3];
    velocity[axis] = if inlet_high { -2.0 } else { 2.0 };
    velocity[(axis + 1) % 3] = 0.25;
    let mut exterior = velocity;
    exterior[axis] = 0.0;
    NsfGrid3::new(
        NsfSpec {
            dims,
            periodic,
            gas,
            temperature_k: 300.0,
            pressure_pa: 100_000.0,
            dx_m: 5e-5,
            temperature_window_k: [250.0, 600.0],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "manufactured thermal-contact fixture",
            maximum_cells: 72,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |_| (300.0 / 300.4, velocity.map(|v| v / unit), 300.4),
    )
    .unwrap()
    .with_characteristic_axis(NsfCharacteristicAxis {
        axis,
        maximum_relative_perturbation: 0.05,
        faces: [NsfCharacteristicBoundary {
            reference_pressure_pa: 100_000.0,
            reference_temperature_k: 300.0,
            reference_velocity_m_s: exterior,
            incoming_characteristic_m_s: 0.0,
        }; 2],
    })
    .unwrap()
}

fn thermal(p: f64, t: f64) -> NsfBoundaryDrive {
    NsfBoundaryDrive::PressureTemperature {
        pressure_offset_pa: p,
        inflow_temperature_k: t,
    }
}

#[test]
fn g1_enthalpy_contact_drives_original_face_energy_flux_on_every_axis_and_side() {
    for axis in 0..3 {
        for high in [false, true] {
            let patch = NsfCharacteristicPatch {
                start: [1, 1],
                end: [2, 2],
            };
            let mut grid = grid(axis, high)
                .with_characteristic_thermal_policy(
                    [NsfCharacteristicThermalPolicy::EnthalpyFlux; 2],
                )
                .unwrap()
                .with_characteristic_patch(high, patch, thermal(8.0, 303.0))
                .unwrap()
                .with_characteristic_patch(!high, patch, thermal(-3.0, 308.0))
                .unwrap();
            let mut controls = [thermal(-3.0, 309.0); 2];
            controls[usize::from(high)] = thermal(8.0, 304.0);
            grid.set_boundary_drives(controls).unwrap();
            let (rho, _) = grid.macro_fields(&grid.f);
            let stress = vec![[0.0; 6]; grid.n];
            for face_high in [false, true] {
                let direction = if face_high { 1 } else { -1 };
                for trace in grid.open_boundary_traces(face_high).unwrap() {
                    let i = grid.index(trace.cell);
                    let face = grid
                        .characteristic_face(grid.open_faces.as_ref(), i, axis, direction)
                        .unwrap();
                    let normal = direction as f64 * face.u[axis];
                    let (flux, _, _) = grid
                        .open_energy_flux(face, i, axis, direction, &rho, &grid.theta, &stress, &[])
                        .unwrap();
                    let kinetic = 0.5 * face.u.iter().map(|v| v * v).sum::<f64>();
                    let expected_h = if face_high == high {
                        assert!(normal < 0.0);
                        let upstream = if patch.contains(trace.cell, axis) {
                            303.0
                        } else {
                            304.0
                        };
                        assert!(trace.temperature_k > 300.4 && trace.temperature_k < upstream);
                        let theta = grid.law.theta_of_kelvin(upstream);
                        theta / (grid.law.gamma - 1.0) + theta
                    } else {
                        assert!(normal > 0.0);
                        assert_eq!(face.theta.to_bits(), grid.theta[i].to_bits());
                        grid.theta[i] / (grid.law.gamma - 1.0) + grid.theta[i]
                    };
                    let expected = face.rho * normal * (expected_h + kinetic);
                    assert!(
                        (direction as f64 * flux - expected).abs()
                            < 2e-13 * expected.abs().max(1e-10)
                    );
                }
            }
        }
    }
}

#[test]
fn g1_constant_material_quadratic_and_zero_flow_limit_are_independent_of_the_root() {
    let mut grid = grid(0, false);
    let adjacent = grid.law.theta_of_kelvin(300.4);
    let k = grid.law.kappa_lattice(adjacent);
    // Freeze the original table only in this manufactured reference.
    // Production always queries the retained material/table.
    grid.law.kappa.fill(k);
    let pressure = 100_000.0;
    let p_lu = pressure / (grid.law.rho_ref_kg_m3 * grid.law.specific_energy_unit_j_kg());
    let cp_over_r = grid.law.gamma / (grid.law.gamma - 1.0);
    // Heating, cooling and equal-temperature cases exercise both bracket orders.
    for upstream_k in [296.0, 300.4, 304.0] {
        let upstream = grid.law.theta_of_kelvin(upstream_k);
        for normal in [-1e-5, -1e-3, -1e-2] {
            let a = p_lu * -normal * cp_over_r;
            let g = 2.0 * k;
            let b = a - g * adjacent;
            // G theta_f^2 + (A-G theta_i) theta_f - A theta_up = 0.
            let reference = (-b + (b * b + 4.0 * g * a * upstream).sqrt()) / (2.0 * g);
            let actual = grid
                .enthalpy_flux_theta(adjacent, pressure, pressure, normal, upstream_k)
                .unwrap();
            assert!((actual - reference).abs() < 3e-15);
        }
    }
    for normal in [0.0, -0.0, 1e-12, 0.01] {
        assert_eq!(
            grid.enthalpy_flux_theta(adjacent, pressure, pressure, normal, 304.0)
                .unwrap()
                .to_bits(),
            adjacent.to_bits()
        );
    }
    let mut previous = f64::INFINITY;
    for normal in [-1e-3, -1e-6, -1e-9, -1e-12] {
        let distance = (grid
            .enthalpy_flux_theta(adjacent, pressure, pressure, normal, 304.0)
            .unwrap()
            - adjacent)
            .abs();
        assert!(distance <= previous);
        previous = distance;
    }
    assert!(previous < 1e-10);
}

#[test]
fn g0_g4_g5_policy_keeps_legacy_bits_and_original_trial_rollback_and_sgs_refusals() {
    let mut original = grid(0, false);
    let mut explicit = grid(0, false)
        .with_characteristic_thermal_policy([NsfCharacteristicThermalPolicy::UpwindEntropy; 2])
        .unwrap();
    for _ in 0..2 {
        assert_eq!(original.state_hash(), explicit.state_hash());
        assert_eq!(original.step().unwrap(), explicit.step().unwrap());
    }
    assert!(
        original
            .with_characteristic_thermal_policy([NsfCharacteristicThermalPolicy::EnthalpyFlux; 2])
            .is_err()
    );
    let policies = [
        NsfCharacteristicThermalPolicy::EnthalpyFlux,
        NsfCharacteristicThermalPolicy::UpwindEntropy,
    ];
    let mut selected = grid(0, false)
        .with_characteristic_thermal_policy(policies)
        .unwrap();
    let before = (
        selected.state_hash(),
        selected.boundary_drives(),
        selected.total_energy_j(),
    );
    let controls = [thermal(2.0, 300.5), NsfBoundaryDrive::PressureOffset(0.0)];
    let work = selected.cell_count();
    let candidate = selected
        .preview_boundary_interval(controls, 1, work, || false)
        .unwrap();
    assert_eq!(candidate.grid().characteristic_thermal_policy(), policies);
    let expected = candidate.grid().state_hash();
    let energy_out: f64 = candidate
        .boundary_exchange()
        .iter()
        .map(|f| f.total_energy_out_j)
        .sum();
    assert!(
        (candidate.grid().total_energy_j() - candidate.previous_energy_j() + energy_out).abs()
            < 5e-12 * candidate.previous_energy_j()
    );
    drop(candidate);
    assert_eq!(
        (
            selected.state_hash(),
            selected.boundary_drives(),
            selected.total_energy_j()
        ),
        before
    );
    assert!(
        selected
            .preview_boundary_interval(controls, 1, work, || true)
            .is_err()
    );
    assert_eq!(
        (
            selected.state_hash(),
            selected.boundary_drives(),
            selected.total_energy_j()
        ),
        before
    );
    selected
        .preview_boundary_interval(controls, 1, work, || false)
        .unwrap()
        .accept();
    assert_eq!(selected.state_hash(), expected);
    assert_eq!(selected.characteristic_thermal_policy(), policies);
    let config = |cs| NsfSgsConfig {
        closure: SgsConfig {
            model: SgsModel::Smagorinsky { cs },
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.9,
    };
    assert!(
        grid(0, false)
            .with_sgs(config(0.17))
            .unwrap()
            .with_characteristic_thermal_policy(policies)
            .is_err()
    );
    assert!(
        grid(0, false)
            .with_characteristic_thermal_policy(policies)
            .unwrap()
            .with_sgs(config(0.17))
            .is_err()
    );
    let molecular = grid(0, false)
        .with_characteristic_thermal_policy(policies)
        .unwrap();
    let hash = molecular.state_hash();
    let zero = molecular.with_sgs(config(0.0)).unwrap();
    assert_eq!(zero.state_hash(), hash);
    assert_eq!(zero.characteristic_thermal_policy(), policies);
}
