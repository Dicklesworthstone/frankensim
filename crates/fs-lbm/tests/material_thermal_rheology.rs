//! G1/G3: the shared SI thermal law drives native LBM collision and heat.
//! Parameters are declared synthetic calibration data, not measured materials.

use fs_exec::CancelGate;
use fs_lbm::material::transport::rheology::RheologySolve;
use fs_lbm::material::{MaterialSlab3, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::rheology::Rheology;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const SHEAR_CELL: usize = 6;
const STRESS_TAU: f64 = 0.8;
const LATTICE_RATE: f64 = 0.02;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1000.0, mu).unwrap()
}
fn controls() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 100,
        relaxation: 0.75,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-12,
    }
}
fn limits() -> HeatSourceLimits {
    HeatSourceLimits {
        max_cells: 16,
        temperature_range_k: [290.0, 330.0],
        balance_tolerance_j: 1e-8,
    }
}
fn law(activation: f64, n: f64) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(1.0, 0.2, 0.5, n, 2.0, [0.0, 1000.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: activation,
            },
            300.0,
            [290.0, 330.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}
fn slab() -> MaterialSlab3 {
    let mut slab = plan_scaling_from_material(&fluid(1.0), DX, DT, 0.1)
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
    // Independent second-order Hermite projection at rho=1, no body force.
    // D_xy=D_yx=gamma/2; its non-equilibrium stress is -2 rho cs^2 tau D.
    // This prescribes gamma*tau, independently of the nonlinear solver.
    let equilibrium = fs_lbm::d3q19::equilibrium3(1.0, [0.03, -0.02, 0.01]);
    slab.lattice.grid.f[SHEAR_CELL] = std::array::from_fn(|q| {
        let e = fs_lbm::d3q19::E3[q];
        equilibrium[q]
            - 3.0
                * fs_lbm::d3q19::W3[q]
                * STRESS_TAU
                * LATTICE_RATE
                * f64::from(e.0)
                * f64::from(e.1)
    });
    // The previous tau is only an initial iterate, not the constitutive answer.
    slab.lattice.grid.tau[SHEAR_CELL] = 1.7;
    slab
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn evolve(model: ThermalCarreauYasuda, original_mu: f64) -> MaterialSlab3 {
    let mut slab = slab();
    let step = slab
        .step_with_material_rheology_and_viscous_heating(
            &[0.0; 16],
            &[[0.0; 3]; 16],
            limits(),
            controls(),
            &CancelGate::new(),
            |point, _| {
                model
                    .evaluate(point.temperature_k, point.shear_rate_per_s)?
                    .bind(fluid(original_mu))
            },
        )
        .unwrap();
    for sample in step
        .transport
        .applied
        .iter()
        .chain(&step.transport.endpoint)
    {
        let retained = sample.material.source();
        let rheology = retained.rheology();
        assert_eq!(rheology.law(), &model);
        assert_eq!(rheology.temperature().temperature_k(), sample.temperature_k);
        assert_eq!(
            rheology.shear_rate_per_s(),
            sample.material.point().shear_rate_per_s
        );
        assert_eq!(
            retained.original_transport().dynamic_viscosity_pa_s(),
            original_mu
        );
        assert_eq!(
            sample.transport.dynamic_viscosity_pa_s(),
            rheology.dynamic_viscosity_pa_s()
        );
        assert_eq!(sample.transport.density_kg_m3(), 1000.0);
        assert_eq!(sample.transport.specific_heat_j_kg_k(), 1000.0);
        assert_eq!(sample.transport.thermal_conductivity_w_m_k(), 1000.0);
        near(
            sample.momentum_relaxation,
            0.5 + 3.0 * rheology.dynamic_viscosity_pa_s() / 1000.0 * DT / (DX * DX),
            2e-12,
        );
    }
    let applied = step
        .transport
        .applied
        .iter()
        .find(|s| s.cell == SHEAR_CELL)
        .unwrap();
    let rate = applied.material.point().shear_rate_per_s;
    near(
        rate,
        STRESS_TAU * LATTICE_RATE / (DT * applied.momentum_relaxation),
        2e-11,
    );
    near(
        step.viscous.power_w_m3[SHEAR_CELL],
        applied.transport.dynamic_viscosity_pa_s() * rate * rate,
        2e-12,
    );
    assert!(step.viscous.power_w_m3[SHEAR_CELL] > 0.0);
    slab
}

#[test]
fn g3_temperature_and_shear_calibration_change_actual_momentum_and_heat() {
    let weak_shift = evolve(law(1000.0, 1.0), 1.0);
    let strong_shift = evolve(law(20_000.0, 1.0), 1.0);
    let shear_thinning = evolve(law(20_000.0, 0.5), 1.0);
    assert_ne!(weak_shift.lattice.grid.f, strong_shift.lattice.grid.f);
    assert_ne!(strong_shift.lattice.grid.f, shear_thinning.lattice.grid.f);
    assert_ne!(strong_shift.lattice.gpop, shear_thinning.lattice.gpop);

    // The declared constitutive viscosity replaces the original mu exactly
    // once; changing only that retained audit value cannot scale it again.
    let different_original = evolve(law(20_000.0, 0.5), 100.0);
    assert_eq!(
        shear_thinning.lattice.grid.f,
        different_original.lattice.grid.f
    );
    assert_eq!(shear_thinning.lattice.gpop, different_original.lattice.gpop);
}

#[test]
fn g1_shared_carreau_special_case_matches_existing_lattice_units() {
    // The legacy lattice owner also accepts thickening and has its own tau
    // floor/cap policy. Compare the common a=2 shear-thinning domain after
    // explicit SI-to-lattice conversion, without assigning SI units to it.
    let curve = CarreauYasuda::new(1.0, 0.2, 0.5, 0.5, 2.0, [0.0, 1e8]).unwrap();
    for (dx, dt, rho) in [(0.01, 0.01, 1000.0), (0.004, 0.0001, 850.0)] {
        let viscosity_scale = dt / (rho * dx * dx);
        let lattice = Rheology::Carreau {
            nu0: viscosity_scale,
            nu_inf: 0.2 * viscosity_scale,
            lambda: 0.5 / dt,
            n: 0.5,
        };
        for rate in [0.0, 0.01, 1.0, 100.0, 1e8] {
            let physical = curve.viscosity_pa_s(rate).unwrap();
            near(
                lattice.viscosity(rate * dt) / viscosity_scale,
                physical,
                3e-13,
            );
        }
    }
}

#[test]
fn g3_calibration_refusal_keeps_native_flow_and_heat_state_unchanged() {
    let mut slab = slab();
    let before_f = slab.lattice.grid.f.clone();
    let before_g = slab.lattice.gpop.clone();
    let before_tau = slab.lattice.grid.tau.clone();
    let model = ThermalCarreauYasuda::new(
        *law(20_000.0, 0.5).curve(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 20_000.0,
            },
            300.0,
            [300.0, 301.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let refusal = slab.step_with_material_rheology_and_viscous_heating(
        &[0.0; 16],
        &[[0.0; 3]; 16],
        limits(),
        controls(),
        &CancelGate::new(),
        |point, _| {
            model
                .evaluate(point.temperature_k, point.shear_rate_per_s)?
                .bind(fluid(1.0))
        },
    );
    assert!(refusal.is_err());
    assert_eq!(slab.lattice.grid.f, before_f);
    assert_eq!(slab.lattice.gpop, before_g);
    assert_eq!(slab.lattice.grid.tau, before_tau);
}
