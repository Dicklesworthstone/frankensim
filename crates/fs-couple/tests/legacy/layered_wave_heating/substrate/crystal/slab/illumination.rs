//! Independent finite-film factors deposit expected heat exactly once.
use super::*;
use fs_couple::layered_wave_heating::substrate::crystal::PreparedCrystalSlabIlluminationHeating;
use fs_em::plane_wave::uniaxial::oriented::absorbing::slab::{
    AbsorbingSlabCoherentIllumination as Group, AbsorbingSlabIllumination as Bank,
};

#[test]
fn unpolarized_multifrequency_film_heat_drives_native_thermal_step() {
    let incident = state(&card("ensemble film incidence", 1.0, 0.0), false, 100_000.0);
    let exit = state(&card("ensemble film exit", 4.0, 0.0), false, 100_000.0);
    let material = crystal_card("ensemble heated film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&exit),
        2.0 * WAVELENGTH_M,
    )
    .unwrap();
    let omega = core::f64::consts::TAU * C0 / WAVELENGTH_M;
    let second_omega = core::f64::consts::TAU * C0 / 650e-9;
    let a = 1e5 * core::f64::consts::FRAC_1_SQRT_2;
    let b = 0.6e5 * core::f64::consts::FRAC_1_SQRT_2;
    // Equal transverse p/s powers at each line require tangential Ex=c*A,
    // Ey=A. Both the polarization weight and line brightness are explicit.
    let fields = [
        film.prepare_absorbing_field(omega, 0.6, [C64::from_re(0.6 * a), C64::ZERO])
            .unwrap(),
        film.prepare_absorbing_field(omega, 0.6, [C64::ZERO, C64::from_re(a)])
            .unwrap(),
        film.prepare_absorbing_field(second_omega, 0.8, [C64::from_re(0.8 * b), C64::ZERO])
            .unwrap(),
        film.prepare_absorbing_field(second_omega, 0.8, [C64::ZERO, C64::from_re(b)])
            .unwrap(),
    ];
    let groups: [_; 4] =
        core::array::from_fn(|i| Group::try_new(&film, &fields[i..i + 1], 1, || false).unwrap());
    let illumination = Bank::try_new(&film, &groups, 4, 4, || false).unwrap();
    let width = 1.5 * WAVELENGTH_M;
    let native = body(&[&material], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let heat = PreparedCrystalSlabIlluminationHeating::try_new(
        &bound,
        &illumination,
        0.25 * WAVELENGTH_M,
        limits(),
        || false,
    )
    .unwrap();
    assert!(std::ptr::eq(heat.source(), &illumination));
    assert!(std::ptr::eq(heat.thermal_body(), &bound));
    assert_eq!(heat.preparation_mode_cell_pairs(), 16 * 4 * (4 + 2));
    let mut powers = Vec::new();
    for (i, cell) in heat.cells().iter().enumerate() {
        assert_eq!(cell.cell_id, bound.cell_ids()[i]);
        // Independent depth quadrature of EACH original Maxwell field.
        let expected = fields
            .iter()
            .map(|field| film_density_integral(field, cell.depth_m[0], cell.depth_m[1]))
            .sum::<f64>();
        assert!(expected > 0.0);
        near(cell.absorbed_w, expected, 2e-9 * expected);
        powers.push(expected + authored[i]);
    }
    assert!(heat.power().residual_w.abs() < 4e-11 * heat.power().entering_w);
    let dt = 2e-9;
    let step = heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(4).unwrap(),
            || false,
        )
        .unwrap();
    let temperatures = thermal_reference(&powers, width, AREA, dt);
    for (i, cell) in heat.cells().iter().enumerate() {
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            temperatures[i],
            3e-8,
        );
        near(step.cell_energy_j[i].1, cell.absorbed_w * dt, 1e-14);
    }
    near(
        step.heat_report.energy.external_heat_j,
        step.optical.absorbed_j + authored.iter().sum::<f64>() * dt,
        1e-9,
    );
    near(
        step.heat_report.energy.stored_energy_change_j,
        step.heat_report.energy.external_heat_j,
        2e-9,
    );
    assert!(matches!(
        heat.advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(4).unwrap(),
            || true,
        ),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}

#[test]
fn film_ensemble_heat_refuses_coherent_beam_groups_and_preserves_dark_evolution() {
    let incident = state(&card("dark ensemble incidence", 1.0, 0.0), false, 100_000.0);
    let material = crystal_card("dark ensemble film", [0.12, 0.42]);
    let [o, e] = axes(&material, 100_000.0);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&incident),
        UniaxialOpticalResponse::try_new(&o, &e, 410e-9).unwrap(),
        [1.0, 2.0, -0.6],
    )
    .unwrap();
    let film = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&incident),
        2.0 * WAVELENGTH_M,
    )
    .unwrap();
    let omega = core::f64::consts::TAU * C0 / WAVELENGTH_M;
    let fields = [
        film.prepare_absorbing_field(omega, 0.6, [C64::ONE, C64::ZERO])
            .unwrap(),
        film.prepare_absorbing_field_at_tangential_index_in_plane(
            omega,
            -0.3,
            [3.0, 4.0],
            [C64::new(0.3, 0.7), C64::ONE],
        )
        .unwrap(),
    ];
    let groups = [Group::try_new(&film, &fields, 2, || false).unwrap()];
    let coherent = Bank::try_new(&film, &groups, 1, 2, || false).unwrap();
    let width = WAVELENGTH_M;
    let native = body(&[&material], &[width], &[4], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    assert!(matches!(PreparedCrystalSlabIlluminationHeating::try_new(
        &bound, &coherent, 0.0, limits(), || false,
    ), Err(LayeredWaveHeatingError::Invalid(message)) if message.contains("one plane wave")));

    let empty_groups = [Group::try_new(&film, &[], 0, || false).unwrap()];
    let dark = Bank::try_new(&film, &empty_groups, 1, 0, || false).unwrap();
    let limits = SubstrateWaveHeatingLimits {
        maximum_mode_cell_pairs: 0,
        ..limits()
    };
    let heat =
        PreparedCrystalSlabIlluminationHeating::try_new(&bound, &dark, 0.0, limits, || false)
            .unwrap();
    assert_eq!(heat.preparation_mode_cell_pairs(), 0);
    assert!(heat.cells().iter().all(|cell| cell.absorbed_w == 0.0));
    assert_eq!(heat.power().entering_w, 0.0);
    assert_eq!(heat.power().leaving_w, 0.0);
    assert_eq!(heat.power().absorbed_w, 0.0);
    assert_eq!(heat.power().residual_w, 0.0);
    // Darkness cannot authorize a thermal column outside its original film.
    assert!(
        PreparedCrystalSlabIlluminationHeating::try_new(
            &bound,
            &dark,
            1.25 * WAVELENGTH_M,
            limits,
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedCrystalSlabIlluminationHeating::try_new(&bound, &dark, 0.0, limits, || true,),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let dt = 2e-9;
    let step = heat
        .advance_to(
            &initial,
            dt,
            conduction(),
            LayeredThermalSolveControls::default(),
            &mut LayeredThermalWorkspace::new(4).unwrap(),
            || false,
        )
        .unwrap();
    let temperatures = thermal_reference(&authored, width, AREA, dt);
    for (i, cell) in heat.cells().iter().enumerate() {
        assert_eq!(step.cell_energy_j[i].1, 0.0);
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            temperatures[i],
            3e-8,
        );
    }
    assert_eq!(step.optical.absorbed_j, 0.0);
    near(
        step.heat_report.energy.external_heat_j,
        authored.iter().sum::<f64>() * dt,
        1e-14,
    );
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
}
