//! Calibrated material kinetics entering the original moving/thermal owners.

use super::{
    Bodies, ImplicitConductionWorkspace, ImplicitThermalScheme, Observe, config, implicit_config,
    laws, near, routes, stiff_graph, stiff_thermal, thermal,
};
use fs_couple::rigid_damping_heat::{
    RigidDampingThermalError,
    temperature::{
        TemperatureDampingLaw, advance_candidate, advance_implicit_candidate,
        shift::{CalibratedKelvinDamping, KelvinTemperatureShiftError},
    },
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    BushingCoordinateScales, PassiveBushingMatrix,
};
use fs_thermal::enthalpy::CellId;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
use std::collections::BTreeMap;

fn arrhenius(validity: [f64; 2]) -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 6000.0 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
        },
        300.0,
        validity,
    )
    .unwrap()
}

#[test]
fn calibrated_shifts_preserve_full_native_forms_and_reference_values() {
    let mut bodies = Bodies::new([0.7, 1.1]);
    let mode = [1.0, -2.0, 0.0, 4.0, 0.0, 0.0];
    bodies.links[0].damping = PassiveBushingMatrix::try_new(
        core::array::from_fn(|i| core::array::from_fn(|j| mode[i] * mode[j])),
        BushingCoordinateScales {
            translation_m: 1.0,
            rotation_rad: 1.0,
        },
    )
    .unwrap();
    let shifts = [
        arrhenius([280.0, 350.0]),
        TemperatureShift::new(
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 8.86,
                c2_k: 101.6,
            },
            310.0,
            [295.0, 350.0],
        )
        .unwrap(),
    ];
    let source = CalibratedKelvinDamping::try_new(&bodies.links, &shifts).unwrap();
    for (edge, reference_temperature) in [300.0, 310.0].into_iter().enumerate() {
        let original = &bodies.links[edge].damping;
        assert_eq!(
            source
                .at_temperature(edge, original, reference_temperature)
                .unwrap(),
            *original
        );
        let temperature = 325.0;
        let expected = if edge == 0 {
            (6000.0_f64 * (1.0 / temperature - 1.0 / 300.0)).exp()
        } else {
            10.0_f64.powf(-8.86 * (temperature - 310.0) / (101.6 + temperature - 310.0))
        };
        let shifted = source.at_temperature(edge, original, temperature).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let actual = shifted.native_coefficients()[i][j];
                let wanted = original.native_coefficients()[i][j] * expected;
                near(actual, wanted, 3e-13 * wanted.abs());
            }
        }
        assert_eq!(shifted.coordinate_scales(), original.coordinate_scales());
    }
    assert_eq!(source.reference_links(), &bodies.links);
    assert_eq!(source.shifts(), &shifts);
}

#[test]
fn calibration_refuses_compounding_an_already_shifted_reference() {
    let bodies = Bodies::new([0.7, 1.1]);
    let shifts = [arrhenius([280.0, 350.0]); 2];
    assert!(CalibratedKelvinDamping::try_new(&bodies.links, &shifts[..1]).is_err());
    let source = CalibratedKelvinDamping::try_new(&bodies.links, &shifts).unwrap();
    let warmed = source
        .at_temperature(0, &bodies.links[0].damping, 320.0)
        .unwrap();
    assert!(matches!(
        source.at_temperature(0, &warmed, 330.0),
        Err(KelvinTemperatureShiftError::Reference { connector: 0 })
    ));
    assert!(matches!(
        source.at_temperature(2, &bodies.links[0].damping, 300.0),
        Err(KelvinTemperatureShiftError::Reference { connector: 2 })
    ));
    assert!(matches!(
        source.at_temperature(0, &bodies.links[0].damping, 351.0),
        Err(KelvinTemperatureShiftError::Shift(_))
    ));
    let original = source
        .at_temperature(0, &bodies.links[0].damping, 300.0)
        .unwrap();
    assert_eq!(original, bodies.links[0].damping);
}

#[test]
fn arrhenius_self_heating_converges_with_original_caloric_storage() {
    let bodies = Bodies::new([1.0, 0.0]);
    let mut calorics = laws();
    calorics[0].a = 0.025; // 2-kg original receiver: total C = 0.05 J/K.
    calorics[0].b = 0.0;
    let initial = thermal(&calorics, false, true);
    let fractions = [BTreeMap::from([(CellId(0), 1.0)]), BTreeMap::new()];
    let shifts = [arrhenius([290.0, 350.0]); 2];
    let source = CalibratedKelvinDamping::try_new(&bodies.links, &shifts).unwrap();

    // Independent physical-coordinate scalar ODE. Reduced mass is 1kg and
    // u0=2m/s. Conservation gives T(u)=300+(2-u^2/2)/0.05, leaving one nonlinear
    // ODE u'=-c(T(u))*u. The source's reduced-time or mechanical solver is not
    // used by this fine RK4 reference.
    let derivative = |u: f64| {
        let temperature = 300.0 + (2.0 - 0.5 * u * u) / 0.05;
        -(6000.0 * (1.0 / temperature - 1.0 / 300.0)).exp() * u
    };
    let mut exact = 2.0;
    let h = 1.0 / 20_000.0;
    for _ in 0..20_000 {
        let k1 = derivative(exact);
        let k2 = derivative(exact + 0.5 * h * k1);
        let k3 = derivative(exact + 0.5 * h * k2);
        let k4 = derivative(exact + h * k3);
        exact += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    let mut errors = Vec::new();
    for intervals in [8, 16, 32] {
        let mut controls = config();
        controls.thermal.max_dt_s = 1.0 / intervals as f64;
        let out = advance_candidate(
            bodies.input(&fractions),
            &initial,
            &(),
            &source,
            &[CellId(0), CellId(0)],
            1.0,
            controls,
            || false,
        )
        .unwrap();
        let slip = out.states()[0].linear_momentum_world().x;
        errors.push((slip - exact).abs());
        let heat = 2.0 - 0.5 * slip * slip;
        near(out.report().dissipated_energy_j, heat, 3e-12);
        near(
            out.thermal().temperature_k(CellId(0)).unwrap(),
            300.0 + heat / 0.05,
            5e-11,
        );
        near(out.report().combined_energy_residual_j, 0.0, 3e-12);
        // Heating reduces the original calibrated viscosity, so more motion
        // remains than under a coefficient frozen at its cold reference.
        assert!(slip > 2.0 * (-1.0_f64).exp() + 0.3);
        let temperature = out.temperatures_k()[0];
        let expected = source
            .at_temperature(0, &bodies.links[0].damping, temperature)
            .unwrap();
        assert_eq!(out.endpoint_links()[0].damping, expected);
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.7 && pair[0] / pair[1] < 2.3,
            "{errors:?}"
        );
    }
    assert!(errors[2] < 0.01, "{errors:?}");
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn endpoint_outside_calibrated_temperature_band_refuses_atomically() {
    let bodies = Bodies::new([0.7, 1.1]);
    let calorics = laws();
    let initial = thermal(&calorics, false, true);
    let fractions = routes();
    let shifts = [arrhenius([300.0, 300.01]); 2];
    let source = CalibratedKelvinDamping::try_new(&bodies.links, &shifts).unwrap();
    let mut controls = config();
    controls.thermal.max_dt_s = 0.5;
    let observer = Observe::default();
    let saved = initial.cells().clone();
    assert!(matches!(
        advance_candidate(
            bodies.input(&fractions),
            &initial,
            &observer,
            &source,
            &[CellId(0), CellId(1)],
            0.5,
            controls,
            || false
        ),
        Err(RigidDampingThermalError::Source {
            error: KelvinTemperatureShiftError::Shift(_),
            ..
        })
    ));
    assert_eq!(initial.cells(), &saved);
    assert_eq!(initial.time_s(), 0.0);
    assert!(observer.intervals.is_empty());
}

#[test]
fn implicit_stiff_heat_calibrates_kelvin_endpoints_and_continues_original_references() {
    let bodies = Bodies::new([0.7, 1.1]);
    let calorics = laws();
    let initial = stiff_thermal(&calorics, 8);
    let graph = stiff_graph(&calorics, &initial);
    let fractions = routes();
    let input = bodies.input(&fractions);
    let shifts = [arrhenius([280.0, 350.0]); 2];
    let source = CalibratedKelvinDamping::try_new(&bodies.links, &shifts).unwrap();
    let cells = [CellId(0), CellId(1)];
    let observer = Observe::default();
    for scheme in [
        ImplicitThermalScheme::BackwardEuler,
        ImplicitThermalScheme::Sdirk2,
    ] {
        let controls = implicit_config(scheme);
        let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
        let (whole, native) = advance_implicit_candidate(
            input,
            &initial,
            &observer,
            &source,
            &cells,
            0.5,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(native.energy.substeps, 2);
        assert_eq!(native.energy, whole.report().thermal);
        for (i, cell) in cells.into_iter().enumerate() {
            let temperature = whole.thermal().temperature_k(cell).unwrap();
            assert!(temperature > 300.0);
            assert_eq!(whole.temperatures_k()[i], temperature);
            let factor = (6000.0 * (1.0 / temperature - 1.0 / 300.0)).exp();
            near(
                whole.endpoint_links()[i].damping.native_coefficients()[0][0],
                [0.7, 1.1][i] * factor,
                2e-14,
            );
        }
        near(whole.report().combined_energy_residual_j, 0.0, 1e-8);

        let (first, _) = advance_implicit_candidate(
            input,
            &initial,
            &observer,
            &source,
            &cells,
            0.25,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        let continued = fs_couple::rigid_damping_heat::RigidDampingThermalInput {
            time_s: first.report().time_s,
            states: first.states(),
            ..input
        };
        // Keep the same original reference links and source calibration. The
        // endpoint forms are diagnostics, not references to shift a second time.
        let (second, _) = advance_implicit_candidate(
            continued,
            first.thermal(),
            first.observer(),
            &source,
            &cells,
            0.5,
            &graph,
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
        assert_eq!(second.states(), whole.states());
        assert_eq!(second.thermal().cells(), whole.thermal().cells());
        assert_eq!(second.endpoint_links(), whole.endpoint_links());
        near(
            first.report().dissipated_energy_j + second.report().dissipated_energy_j,
            whole.report().dissipated_energy_j,
            1e-12,
        );

        let narrow = [arrhenius([300.0, 300.001]); 2];
        let narrow_source = CalibratedKelvinDamping::try_new(&bodies.links, &narrow).unwrap();
        assert!(matches!(
            advance_implicit_candidate(
                input,
                &initial,
                &observer,
                &narrow_source,
                &cells,
                0.25,
                &graph,
                controls,
                &mut workspace,
                || false
            ),
            Err(RigidDampingThermalError::Source {
                error: KelvinTemperatureShiftError::Shift(_),
                ..
            })
        ));
    }
    assert_eq!(source.reference_links(), &bodies.links);
    assert_eq!(initial.time_s(), 0.0);
    assert!(
        initial
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    assert_eq!(observer, Observe::default());
}
