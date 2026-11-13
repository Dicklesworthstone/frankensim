//! G1/G3/G4: chemical and radiant heat reach the original native gas traction.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceInterfaceFace, FreeSurfaceInterfaceSolve,
};

// The rigid 8^3 vessel has 4*6*6 liquid cells of side 1 mm. This physical
// atmosphere volume is fixed; numerical density/fill never supplies a new one.
const GAS_VOLUME: f64 = 3.68e-7;

fn parcel(law: &GasReaction) -> IsochoricGasParcel {
    let initial = gas(law);
    IsochoricGasParcel::new(initial.model().clone(), initial.point().clone(), GAS_VOLUME).unwrap()
}

fn front_velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> f64 {
    liquid
        .flow()
        .sample_at_si(3, 3, 3)
        .unwrap()
        .unwrap()
        .velocity_m_s[0]
}

#[test]
fn g1_chemical_heat_loads_native_momentum_with_post_thermal_pressure() {
    for release in [-2000.0, 2000.0] {
        with_liquid_pressure(false, Some(100000.0), |liquid, gate| {
            let law = reaction(5.0, release, 8.0);
            let original = parcel(&law);
            let next = law
                .pressurized_liquid_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &[face()],
                    0.05,
                    controls(),
                    gate,
                )
                .unwrap();
            let unloaded = law
                .liquid_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &[face()],
                    controls(),
                    gate,
                )
                .unwrap();

            // Independent isomer extent, cv and four Fourier film increments.
            // Formation-energy release changes T; it is not an extra heat debit.
            let extent = DT * GAS_VOLUME * 5.0;
            let amount = 100000.0 * GAS_VOLUME / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0);
            let capacity = 2.5 * 100000.0 * GAS_VOLUME / 400.0;
            let conductivity = 0.03 * (1.0 + 7.0 * extent / amount);
            let conductance = 1e-6 / (0.0005 / 0.5 + 0.0005 / conductivity);
            let intermediate = 400.0 + release * extent / capacity;
            let (mut tl, mut tg) = (400.0, intermediate);
            for _ in 0..4 {
                let heat = DT / 4.0 * conductance * (tg - tl);
                tl += heat / 0.002;
                tg -= heat / capacity;
            }
            let expected_pressure = 100000.0 * tg / 400.0;
            let finite = next
                .report()
                .thermal
                .conduction
                .finite_gas
                .as_ref()
                .unwrap();
            let pressure = finite.pressure.unwrap();
            near(
                pressure.initial_pressure_pa,
                100000.0 * intermediate / 400.0,
                1e-6,
            );
            near(pressure.applied_pressure_pa, expected_pressure, 1e-6);
            assert_eq!(
                pressure.applied_pressure_pa,
                next.gas().state().state().pressure
            );
            assert_eq!(pressure.volume_m3, GAS_VOLUME);
            assert_eq!(pressure.maximum_density_excursion, 0.05);
            assert_ne!(pressure.applied_pressure_pa, pressure.initial_pressure_pa);
            near(next.gas().state().state().temperature, tg, 2e-9);
            assert_eq!(next.gas(), unloaded.gas());

            // Incoming gas links at this flat front have total weight 1/6.
            // Their reconstructed density and momentum give u_x=-d/(1+d),
            // where d=delta_p/(rho_liquid*(dx/dt)^2), and dx/dt=1 m/s.
            let perturbation = (expected_pressure - 100000.0) / 1000.0;
            near(
                front_velocity(next.liquid()),
                -perturbation / (1.0 + perturbation),
                2e-9,
            );
            near(front_velocity(unloaded.liquid()), 0.0, 2e-13);
            assert!(release * front_velocity(next.liquid()) < 0.0);
            assert_ne!(
                next.liquid().flow().lattice().grid().state_bits(),
                unloaded.liquid().flow().lattice().grid().state_bits()
            );
            assert_eq!(next.liquid().accepted_steps(), 1);
            assert_eq!(next.liquid().physical_time_s(), DT);
            near(next.report().supplied_thermal_energy_j, 0.0, 0.0);
            assert!(
                next.report().combined_energy_residual_j.abs()
                    < controls().total_energy_tolerance_j
            );
            assert_eq!(original, parcel(&law));
            assert_eq!(liquid.accepted_steps(), 0);
        });
    }
}

#[test]
fn g3_radiation_and_chemistry_share_one_gas_pressure_and_one_flow_step() {
    with_liquid_pressure(false, Some(100000.0), |liquid, gate| {
        let law = reaction(5.0, 2000.0, 8.0);
        let original = parcel(&law);
        let surface = FreeSurfaceInterfaceFace {
            face: face(),
            emissivity: 0.8.into(),
            surroundings_k: 450.0,
        };
        let next = law
            .pressurized_liquid_surface_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[surface],
                FreeSurfaceInterfaceSolve::default(),
                0.05,
                controls(),
                gate,
            )
            .unwrap();
        let unloaded = law
            .liquid_surface_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[surface],
                FreeSurfaceInterfaceSolve::default(),
                controls(),
                gate,
            )
            .unwrap();
        let report = next.report();
        let exchange = report.thermal.conduction.surface.as_ref().unwrap();
        let pressure = report
            .thermal
            .conduction
            .finite_gas
            .as_ref()
            .unwrap()
            .pressure
            .unwrap();
        assert!(exchange.heat_from_radiation_j > 0.0);
        near(
            report.supplied_thermal_energy_j,
            exchange.heat_from_radiation_j,
            1e-12,
        );
        assert_eq!(next.gas(), unloaded.gas());
        assert_eq!(
            pressure.applied_pressure_pa,
            next.gas().state().state().pressure
        );
        let perturbation = (pressure.applied_pressure_pa - 100000.0) / 1000.0;
        near(
            front_velocity(next.liquid()),
            -perturbation / (1.0 + perturbation),
            2e-9,
        );
        assert!(front_velocity(next.liquid()) < 0.0);
        near(front_velocity(unloaded.liquid()), 0.0, 2e-13);
        assert!(report.combined_energy_residual_j.abs() < controls().total_energy_tolerance_j);
        assert_eq!(next.liquid().accepted_steps(), 1);
        assert_eq!(next.liquid().physical_time_s(), DT);
    });
}

#[test]
fn g4_pressure_refusal_preserves_chemistry_heat_clock_and_retry() {
    with_liquid_pressure(false, Some(100000.0), |liquid, gate| {
        let law = reaction(5.0, 2000.0, 8.0);
        let original = parcel(&law);
        let bits = liquid.flow().lattice().grid().state_bits();
        let enthalpy = liquid.total_enthalpy_j().unwrap();
        let run = |maximum, gate| {
            law.pressurized_liquid_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[face()],
                maximum,
                controls(),
                gate,
            )
        };
        // This fails after chemistry and finite-capacity heat were evaluated.
        let error = run(1e-8, gate)
            .err()
            .expect("pressure excursion must refuse");
        assert!(format!("{error:?}").contains("InvalidGasPressure"));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            run(0.05, &cancelled),
            Err(ReactiveLiquidGasError::Cancelled)
        ));
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(original, parcel(&law));
        let retry = run(0.05, gate).unwrap();
        let uninterrupted = run(0.05, gate).unwrap();
        assert_eq!(retry.gas(), uninterrupted.gas());
        assert_eq!(
            retry.liquid().flow().lattice().grid().state_bits(),
            uninterrupted.liquid().flow().lattice().grid().state_bits()
        );
        assert_eq!(
            retry.liquid().total_enthalpy_j().unwrap(),
            uninterrupted.liquid().total_enthalpy_j().unwrap()
        );
        assert_eq!(retry.liquid().accepted_steps(), 1);
    });
    with_liquid(false, |liquid, gate| {
        let law = reaction(5.0, 2000.0, 8.0);
        let original = parcel(&law);
        let error = law
            .pressurized_liquid_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[face()],
                0.05,
                controls(),
                gate,
            )
            .err()
            .expect("missing physical pressure reference must refuse");
        assert!(format!("{error:?}").contains("pressure reference"));
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(original, parcel(&law));
    });
}
