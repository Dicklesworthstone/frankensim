//! Independently oriented optical layers share their native thermal cells.
use super::*;
use fs_couple::layered_wave_heating::substrate::crystal::PreparedCrystalStackHeating;
use fs_em::plane_wave::uniaxial::oriented::stack::{
    OrientedUniaxialStack, UniaxialStackLayer, UniaxialStackLimits,
};

#[test]
fn distinct_crystal_layers_deposit_tensor_work_once_on_native_thermal_state() {
    let exterior = state(&card("stack exterior", 1.0, 0.0), false, 100_000.0);
    let cards = [
        crystal_card("first heated layer", [0.12, 0.42]),
        crystal_card("second heated layer", [0.28, 0.05]),
    ];
    let principal = cards.each_ref().map(|card| axes(card, 100_000.0));
    let width = 0.75 * WAVELENGTH_M;
    let orientations = [[1.0, 2.0, -0.6], [-0.7, 0.2, 1.0]];
    let layers = core::array::from_fn::<_, 2, _>(|i| {
        UniaxialStackLayer::try_new(
            UniaxialOpticalResponse::try_new(&principal[i][0], &principal[i][1], 410e-9).unwrap(),
            orientations[i],
            width,
        )
        .unwrap()
    });
    let stack = OrientedUniaxialStack::try_new(
        DielectricResponse::Lossless(&exterior),
        DielectricResponse::Lossless(&exterior),
        &layers,
        UniaxialStackLimits {
            maximum_layers: 2,
            maximum_workspace_bytes: 1 << 20,
        },
    )
    .unwrap();
    let field = stack
        .prepare_absorbing_field_in_plane(
            core::f64::consts::TAU * C0 / WAVELENGTH_M,
            0.8,
            [0.6, 0.8],
            [C64::new(1.4e5, 0.4e5), C64::new(-0.6e5, 1.0e5)],
            || false,
        )
        .unwrap();
    let native = body(&[&cards[0], &cards[1]], &[width; 2], &[2, 2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let authored = [1.0, 2.0, 3.0, 4.0];
    let initial = initial(&bound, &authored);
    let original = initial.cells().clone();
    let heat = PreparedCrystalStackHeating::try_new(&bound, &field, limits(), || false).unwrap();
    assert!(core::ptr::eq(heat.source(), &field));
    assert!(core::ptr::eq(heat.thermal_body(), &bound));
    assert_eq!(heat.preparation_mode_cell_pairs(), 16 * (4 + 2));
    let mut powers = Vec::new();
    for (index, cell) in heat.cells().iter().enumerate() {
        assert_eq!(cell.cell_id, bound.cell_ids()[index]);
        let layer = index / 2;
        let local_start = (index % 2) as f64 * width * 0.5;
        let local_end = local_start + width * 0.5;
        let expected = field
            .layer_absorption(layer, AREA, local_start, local_end)
            .unwrap()
            .absorbed_power_w();
        near(cell.absorbed_w, expected, 2e-12 * expected);
        powers.push(expected + authored[index]);
    }
    assert!(heat.power().residual_w.abs() < 3e-11 * heat.power().entering_w);
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
    let expected = thermal_reference(&powers, 2.0 * width, AREA, dt);
    for (index, cell) in heat.cells().iter().enumerate() {
        near(
            step.thermal.temperature_k(cell.cell_id).unwrap(),
            expected[index],
            3e-8,
        );
        near(step.cell_energy_j[index].1, cell.absorbed_w * dt, 1e-14);
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
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells(), &original);
    assert!(matches!(
        PreparedCrystalStackHeating::try_new(&bound, &field, limits(), || true),
        Err(LayeredWaveHeatingError::Cancelled)
    ));
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
    assert_eq!(initial.cells(), &original);

    // Original cards, fixed query context, exact thicknesses and ordered
    // layer count remain requirements even when the temperatures coincide.
    let foreign = body(&[&cards[1], &cards[0]], &[width; 2], &[2, 2], AREA);
    let foreign_laws = foreign.thermal_laws(false).unwrap();
    let foreign_bound = foreign.bind(&foreign_laws, ids()).unwrap();
    assert!(
        PreparedCrystalStackHeating::try_new(&foreign_bound, &field, limits(), || false).is_err()
    );
    let shortened = body(
        &[&cards[0], &cards[1]],
        &[width, width * 0.9],
        &[2, 2],
        AREA,
    );
    let short_laws = shortened.thermal_laws(false).unwrap();
    let short_bound = shortened.bind(&short_laws, ids()).unwrap();
    assert!(
        PreparedCrystalStackHeating::try_new(&short_bound, &field, limits(), || false).is_err()
    );
    assert!(
        PreparedCrystalStackHeating::try_new(
            &bound,
            &field,
            SubstrateWaveHeatingLimits {
                maximum_mode_cell_pairs: 95,
                ..limits()
            },
            || false
        )
        .is_err()
    );
}
