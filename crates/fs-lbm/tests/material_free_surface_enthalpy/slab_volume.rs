//! Actual native interface motion from original thermal density, then native continuation.
use super::*;
use fs_exec::Cx;
use fs_lbm::conjugate::DensityProjectionControls;
use fs_lbm::d3q19::Face3;
use fs_lbm::material::free_surface::MaterialFreeSurfaceVolumeControls3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceSlabVolumeControls, MaterialFreeSurfaceEnthalpy3,
};
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve};

#[path = "slab_volume/heat.rs"]
mod heat;

#[derive(Clone, Copy)]
struct Case {
    face: Face3,
    hot_density: f64,
    start_hot: bool,
    increment_j_kg: f64,
    datum: f64,
    fill: f64,
    depth_cells: usize,
    terminal_fill: f64,
    speed: f64,
    hole: bool,
    thermal_gradient: bool,
}
impl Default for Case {
    fn default() -> Self {
        Self {
            face: Face3::XMin,
            hot_density: 980.0,
            start_hot: false,
            increment_j_kg: 400_000.0,
            datum: 0.0,
            fill: 1.0,
            depth_cells: 2,
            terminal_fill: 1.0,
            speed: 0.0,
            hole: false,
            thermal_gradient: false,
        }
    }
}

fn controls() -> FreeSurfaceSlabVolumeControls {
    FreeSurfaceSlabVolumeControls {
        volume: MaterialFreeSurfaceVolumeControls3 {
            projection: DensityProjectionControls {
                duration_s: 0.001,
                relative_tolerance: 1e-12,
                maximum_volume_defect_m3_s: 1e-18,
                maximum_energy_defect_j: 1e-18,
            },
            max_cells: 512,
            maximum_initial_volume_defect_m3: 1e-20,
            maximum_source_density_defect_kg_m3: 1e-9,
            maximum_enthalpy_defect_j: 1e-8,
        },
        maximum_face_flux_defect_m3_s: 1e-16,
        maximum_mass_defect_kg: 1e-16,
        maximum_momentum_defect_kg_m_s: 1e-16,
        maximum_geometry_volume_defect_m3: 1e-18,
        maximum_native_centroid_offset_m: 0.0005,
        maximum_kinetic_mapping_change_j: 1e-8,
    }
}

fn with_case(
    case: Case,
    run: impl FnOnce(
        &MaterialFreeSurfaceEnthalpy3<'_, '_>,
        &ResolvedLiquidSurfaceState,
        &LiquidTemperatureLaw<'_>,
        &CancelGate,
        &Cx<'_>,
    ),
) {
    with_case_and_original_law(case, |old, source, law, _, gate, cx| {
        run(old, source, law, gate, cx);
    });
}

fn with_case_and_original_law(
    case: Case,
    run: impl FnOnce(
        &MaterialFreeSurfaceEnthalpy3<'_, '_>,
        &ResolvedLiquidSurfaceState,
        &LiquidTemperatureLaw<'_>,
        &LiquidTemperatureLaw<'_>,
        &CancelGate,
        &Cx<'_>,
    ),
) {
    let card = card(0.1, case.hot_density);
    let chart = ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::with_quantity(
                "specific_heat_capacity",
                QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
            ),
            density: &PropertyKey::with_quantity(
                "density",
                QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            ),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(500.0),
            selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        case.datum,
    )
    .unwrap();
    let t0 = if case.start_hot { 500.0 } else { 300.0 };
    let h0 = case.datum + if case.start_hot { 400_000.0 } else { 0.0 };
    let initial = chart.curve().state_at_specific_enthalpy(h0).unwrap();
    let proposed = chart
        .curve()
        .advance_specific_energy(initial, case.increment_j_kg)
        .unwrap();
    let t1 = proposed.temperature_k();
    let original_law = LiquidTemperatureLaw::try_new(
        &card,
        &point(t0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let endpoint_law = LiquidTemperatureLaw::try_new(
        &card,
        &point(t1),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original_source = resolve_liquid_surface_state(
        &card,
        &point(t0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    // The endpoint is the surface state of the endpoint law's whole source
    // bundle, as the library's caloric endpoint builds it: the slab advance
    // requires the endpoint to keep every original claim, heat capacity and
    // conductivity included, which a density and viscosity resolution lacks.
    let endpoint_source = ResolvedLiquidSurfaceState::try_from_material(
        endpoint_law.reference_state().material().clone(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
    )
    .unwrap();
    let axis = case.face as usize / 2;
    let upper = case.face as usize % 2 == 1;
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let coordinates = [x, y, z];
                let normal = if upper {
                    7 - coordinates[axis]
                } else {
                    coordinates[axis]
                };
                if normal < case.depth_cells && !(case.hole && coordinates == [0; 3]) {
                    let mut velocity_m_s = [0.0; 3];
                    velocity_m_s[axis] = case.speed;
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates,
                        fill_fraction: case.fill
                            * if normal + 1 == case.depth_cells {
                                case.terminal_fill
                            } else {
                                1.0
                            },
                        velocity_m_s,
                    });
                }
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&original_source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let thermal = flow
            .with_enthalpy(
                &original_law,
                chart.curve(),
                FreeSurfaceEnthalpyControls::default(),
                |p| {
                    h0 + if case.thermal_gradient && p == [0; 3] {
                        1000.0
                    } else {
                        0.0
                    }
                },
                &gate,
            )
            .unwrap()
            .with_time_origin(0.25)
            .unwrap();
        run(
            &thermal,
            &endpoint_source,
            &endpoint_law,
            &original_law,
            &gate,
            cx,
        );
    });
}

#[test]
fn g3_source_density_moves_native_slab_on_every_axis_and_preserves_material_inventories() {
    for face in Face3::ALL {
        for hot_density in [980.0, 1020.0] {
            with_case(
                Case {
                    face,
                    hot_density,
                    ..Case::default()
                },
                |old, source, law, gate, cx| {
                    let before_bits = old.flow().lattice().grid().state_bits();
                    let before_h = old.total_enthalpy_j().unwrap().to_bits();
                    let step = old
                        .advance_uniform_slab_volume(
                            face,
                            400_000.0,
                            source,
                            law,
                            controls(),
                            gate,
                            cx,
                        )
                        .unwrap();
                    let rho0 = old.flow().scaling().material().density_kg_m3();
                    let rho1 = source.density_kg_m3();
                    let stretch = rho0 / rho1;
                    let height = 0.002 * stretch;
                    let rate = (rho0 - rho1) / (rho1 * 0.001);
                    let mass = 128.0 * old.flow().scaling().mass_unit_kg();
                    let sign = if face as usize % 2 == 0 { 1.0 } else { -1.0 };
                    let axis = face as usize / 2;
                    let momentum = sign * mass * rate * 0.002 / 2.0;
                    let kinetic = mass * (rate * 0.002).powi(2) / 6.0;
                    near(step.report.final_height_m, height, 2e-13);
                    near(step.report.mass_before_kg, mass, 2e-13);
                    near(step.report.mass_after_kg, mass, 2e-13);
                    near(
                        step.report.enthalpy_after_j - step.report.enthalpy_before_j,
                        mass * 400_000.0,
                        2e-12,
                    );
                    near(step.report.momentum_after_kg_m_s[axis], momentum, 2e-12);
                    near(step.report.pressure_impulse_kg_m_s[axis], momentum, 2e-12);
                    near(step.report.affine_kinetic_energy_j, kinetic, 2e-12);
                    assert!(step.report.native_kinetic_energy_j <= kinetic + 1e-18);
                    assert!(step.report.affine_to_native_kinetic_change_j <= 1e-18);
                    assert!(step.report.face_to_affine_kinetic_change_j <= 1e-18);
                    near(
                        step.report.required_mechanical_work_j,
                        step.projection.projection.report.pressure_source_work_j,
                        1e-13,
                    );
                    assert!(
                        step.report.required_mechanical_work_j
                            > step.report.native_kinetic_energy_j
                    );
                    assert!(step.report.mechanical_energy_defect_j.abs() < 1e-18);
                    assert_eq!(step.liquid.accepted_steps(), old.accepted_steps() + 1);
                    near(step.liquid.physical_time_s(), 0.251, 1e-14);
                    // Source density changes the real interface, not numerical pressure.
                    let cells_deep = (height / 0.001).ceil() as usize;
                    let mut at = [3; 3];
                    at[axis] = if sign > 0.0 {
                        cells_deep - 1
                    } else {
                        8 - cells_deep
                    };
                    let sample = step
                        .liquid
                        .flow()
                        .sample_at_si(at[0], at[1], at[2])
                        .unwrap()
                        .unwrap();
                    near(
                        sample.fill_fraction,
                        height / 0.001 - (cells_deep - 1) as f64,
                        2e-12,
                    );
                    near(sample.density_kg_m3, rho1, 2e-14);
                    assert_eq!(sample.pressure_perturbation_pa, 0.0);
                    near(
                        step.liquid.specific_enthalpy_at(at).unwrap(),
                        400_000.0,
                        1e-14,
                    );
                    let mut beyond = at;
                    beyond[axis] = if sign > 0.0 {
                        cells_deep
                    } else {
                        7 - cells_deep
                    };
                    assert_eq!(
                        step.liquid
                            .flow()
                            .region_mass_kg(|x, y, z| [x, y, z] == beyond,)
                            .unwrap(),
                        0.0,
                    );
                    if let Some(sample) = step
                        .liquid
                        .flow()
                        .sample_at_si(beyond[0], beyond[1], beyond[2])
                        .unwrap()
                    {
                        assert_eq!(sample.fill_fraction, 0.0);
                    }
                    assert_eq!(
                        step.report.newly_wet_cells,
                        if hot_density < 1000.0 { 64 } else { 0 }
                    );
                    assert_eq!(old.flow().lattice().grid().state_bits(), before_bits);
                    assert_eq!(old.total_enthalpy_j().unwrap().to_bits(), before_h);
                    assert_eq!(old.accepted_steps(), 0);
                },
            );
        }
    }
}

#[test]
fn g3_interior_heating_cooling_keeps_original_calorics_and_continues_native_enthalpy() {
    for start_hot in [false, true] {
        for magnitude in [0.1, 12345.6789] {
            let dh = if start_hot { -magnitude } else { magnitude };
            with_case(
                Case {
                    start_hot,
                    increment_j_kg: dh,
                    ..Case::default()
                },
                |old, source, law, gate, cx| {
                    let mut step = old
                        .advance_uniform_slab_volume(
                            Face3::XMin,
                            dh,
                            source,
                            law,
                            controls(),
                            gate,
                            cx,
                        )
                        .unwrap();
                    let initial_mass = old.flow().ledger_mass_kg().unwrap();
                    let initial_enthalpy = old.total_enthalpy_j().unwrap();
                    near(
                        step.liquid.flow().ledger_mass_kg().unwrap(),
                        initial_mass,
                        2e-12,
                    );
                    near(
                        step.liquid.total_enthalpy_j().unwrap() - initial_enthalpy,
                        initial_mass * dh,
                        2e-6,
                    );
                    let agreement = step.liquid.density_agreement().unwrap().unwrap();
                    assert!(agreement.maximum_defect_kg_m3 <= agreement.maximum_allowed_kg_m3);
                    let accepted_enthalpy = step.liquid.total_enthalpy_j().unwrap();
                    let tick = step.liquid.step(gate).unwrap();
                    near(tick.mass_after_kg, initial_mass, 2e-12);
                    near(
                        step.liquid.total_enthalpy_j().unwrap(),
                        accepted_enthalpy,
                        2e-12,
                    );
                    assert!(
                        tick.density_agreement.unwrap().maximum_defect_kg_m3
                            <= agreement.maximum_allowed_kg_m3
                    );
                    near(step.liquid.physical_time_s(), 0.252, 1e-14);
                    // Ordinary constant-reference stepping cannot hide another
                    // genuine density change, even though arithmetic agreement is explicit.
                    let p = [0; 3];
                    let cell_mass = step
                        .liquid
                        .flow()
                        .region_mass_kg(|x, y, z| [x, y, z] == p)
                        .unwrap();
                    let before = step.liquid.flow().lattice().grid().state_bits();
                    let clock = step.liquid.physical_time_s().to_bits();
                    let energy = cell_mass * if start_hot { -1000.0 } else { 1000.0 };
                    assert!(matches!(
                        step.liquid.step_with_heat(
                            &[FreeSurfaceCellHeat {
                                coordinates: p,
                                energy_j: energy
                            }],
                            gate,
                        ),
                        Err(FreeSurfaceEnthalpyError::Invalid(
                            "ordinary native step exceeds its fixed-reference density accuracy budget"
                        ))
                    ));
                    assert_eq!(step.liquid.flow().lattice().grid().state_bits(), before);
                    assert_eq!(step.liquid.physical_time_s().to_bits(), clock);
                },
            );
        }
    }
}

#[test]
fn g3_zero_net_heat_never_manufactures_motion_and_enthalpy_datum_does_not_set_expansion() {
    with_case(
        Case {
            increment_j_kg: 0.0,
            ..Case::default()
        },
        |old, source, law, gate, cx| {
            let step = old
                .advance_uniform_slab_volume(Face3::XMin, 0.0, source, law, controls(), gate, cx)
                .unwrap();
            assert_eq!(
                old.flow().lattice().grid().state_bits(),
                step.liquid.flow().lattice().grid().state_bits()
            );
            assert_eq!(step.report.native_kinetic_energy_j, 0.0);
            assert_eq!(step.report.required_mechanical_work_j, 0.0);
            assert_eq!(step.report.momentum_after_kg_m_s, [0.0; 3]);
            assert_eq!(step.report.final_height_m, 0.002);
            assert_eq!(step.report.mass_before_kg, step.report.mass_after_kg);
        },
    );
    let mut baseline = None;
    for datum in [0.0, 1e8] {
        with_case(
            Case {
                datum,
                ..Case::default()
            },
            |old, source, law, gate, cx| {
                let step = old
                    .advance_uniform_slab_volume(
                        Face3::XMin,
                        400_000.0,
                        source,
                        law,
                        controls(),
                        gate,
                        cx,
                    )
                    .unwrap();
                let bits = step.liquid.flow().lattice().grid().state_bits();
                if let Some(expected) = &baseline {
                    assert_eq!(&bits, expected);
                } else {
                    baseline = Some(bits);
                }
                assert!(
                    step.report.enthalpy_defect_j.abs()
                        <= controls().volume.maximum_enthalpy_defect_j
                );
            },
        );
    }
}

#[test]
fn g0_unsupported_slab_state_source_and_cancellation_preserve_original_prefix() {
    for case in [
        Case {
            fill: 0.5,
            ..Case::default()
        },
        Case {
            speed: 0.001,
            ..Case::default()
        },
        Case {
            depth_cells: 7,
            terminal_fill: 0.95,
            ..Case::default()
        },
        Case {
            hole: true,
            ..Case::default()
        },
        Case {
            hot_density: 1000.0,
            thermal_gradient: true,
            ..Case::default()
        },
    ] {
        with_case(case, |old, source, law, gate, cx| {
            let bits = old.flow().lattice().grid().state_bits();
            let h = old.total_enthalpy_j().unwrap().to_bits();
            assert!(
                old.advance_uniform_slab_volume(
                    case.face,
                    case.increment_j_kg,
                    source,
                    law,
                    controls(),
                    gate,
                    cx
                )
                .is_err()
            );
            assert_eq!(old.flow().lattice().grid().state_bits(), bits);
            assert_eq!(old.total_enthalpy_j().unwrap().to_bits(), h);
            assert_eq!(old.accepted_steps(), 0);
        });
    }
    with_case(Case::default(), |old, source, law, gate, cx| {
        let before = old.flow().lattice().grid().state_bits();
        assert!(
            old.advance_uniform_slab_volume(
                Face3::XMin,
                400_000.0,
                old.flow().scaling().material(),
                law,
                controls(),
                gate,
                cx
            )
            .is_err()
        );
        let strict_mapping = FreeSurfaceSlabVolumeControls {
            maximum_kinetic_mapping_change_j: 0.0,
            ..controls()
        };
        assert!(
            old.advance_uniform_slab_volume(
                Face3::XMin,
                400_000.0,
                source,
                law,
                strict_mapping,
                gate,
                cx,
            )
            .is_err()
        );
        gate.request();
        assert!(
            old.advance_uniform_slab_volume(
                Face3::XMin,
                400_000.0,
                source,
                law,
                controls(),
                gate,
                cx
            )
            .is_err()
        );
        assert_eq!(old.flow().lattice().grid().state_bits(), before);
        assert_eq!(old.accepted_steps(), 0);
    });
}

#[test]
fn g3_large_enthalpy_datum_retains_a_small_thermal_increment_before_subtraction() {
    with_case(
        Case {
            datum: 1e8,
            increment_j_kg: 0.1,
            ..Case::default()
        },
        |old, source, law, gate, cx| {
            let step = old
                .advance_uniform_slab_volume(Face3::XMin, 0.1, source, law, controls(), gate, cx)
                .unwrap();
            let mass = old.flow().ledger_mass_kg().unwrap();
            near(
                step.report.requested_enthalpy_increment_j,
                mass * 0.1,
                2e-14,
            );
            assert!(
                step.report.enthalpy_defect_j.abs() <= controls().volume.maximum_enthalpy_defect_j
            );
            assert!(step.report.mass_defect_enthalpy_datum_j.is_finite());
            assert!(step.report.final_height_m > 0.002);
            assert!(step.report.momentum_after_kg_m_s[0] > 0.0);
            near(step.liquid.physical_time_s(), 0.251, 1e-14);
        },
    );
}

#[test]
fn g3_partial_terminal_and_subcell_slabs_use_the_actual_pressure_geometry() {
    for face in Face3::ALL {
        for (depth_cells, terminal_fill) in [(2, 0.35), (1, 0.6)] {
            for hot_density in [980.0, 1020.0] {
                with_case(
                    Case {
                        face,
                        depth_cells,
                        terminal_fill,
                        hot_density,
                        ..Case::default()
                    },
                    |old, source, law, gate, cx| {
                        let step = old
                            .advance_uniform_slab_volume(
                                face,
                                400_000.0,
                                source,
                                law,
                                controls(),
                                gate,
                                cx,
                            )
                            .unwrap();
                        let axis = face as usize / 2;
                        let upper = face as usize % 2 == 1;
                        let depth = (depth_cells - 1) as f64 + terminal_fill;
                        let mass = 64.0 * depth * old.flow().scaling().mass_unit_kg();
                        let stretch = 1000.0 / hot_density;
                        let height = depth * 0.001;
                        let rate = (1000.0 - hot_density) / hot_density / 0.001;
                        near(step.report.initial_height_m, height, 2e-13);
                        near(step.report.final_height_m, height * stretch, 2e-13);
                        near(step.report.mass_after_kg, mass, 2e-13);
                        near(
                            step.report.affine_kinetic_energy_j,
                            mass * (rate * height).powi(2) / 6.0,
                            2e-12,
                        );
                        let terminal = if upper {
                            8 - depth_cells
                        } else {
                            depth_cells - 1
                        };
                        let mut at = [3; 3];
                        at[axis] = terminal;
                        let c = step.projection.domain.index(at[0], at[1], at[2]);
                        near(
                            step.projection.domain.width(axis, terminal),
                            terminal_fill * 0.001,
                            2e-13,
                        );
                        near(
                            step.projection.domain.volume(c),
                            terminal_fill * 1e-9,
                            2e-13,
                        );
                        let gas_face = 2 * axis + usize::from(!upper);
                        let gas_flux = step.projection.projection.field.outward(
                            &step.projection.domain,
                            c,
                            gas_face,
                        );
                        near(gas_flux, rate * height * 1e-6, 2e-11);
                        assert!(
                            step.report.maximum_native_centroid_offset_m
                                >= 0.5 * (1.0 - terminal_fill) * 0.001
                        );
                        assert!(step.report.mechanical_energy_defect_j.abs() <= 1e-18);
                        assert_eq!(step.report.momentum_before_kg_m_s, [0.0; 3]);
                        assert_eq!(step.report.native_kinetic_energy_before_j, 0.0);
                        let too_small = FreeSurfaceSlabVolumeControls {
                            maximum_native_centroid_offset_m: 0.0,
                            ..controls()
                        };
                        assert!(
                            old.advance_uniform_slab_volume(
                                face, 400_000.0, source, law, too_small, gate, cx,
                            )
                            .is_err()
                        );
                    },
                );
            }
        }
    }
}

#[test]
fn g3_retained_native_motion_repeats_heating_and_cooling_with_actual_pressure_work() {
    for face in [Face3::XMin, Face3::YMax, Face3::ZMin] {
        with_case_and_original_law(
            Case {
                face,
                terminal_fill: 0.7,
                ..Case::default()
            },
            |old, hot_source, hot_law, cold_law, gate, cx| {
                // Back at the original temperature, the endpoint is again the
                // law's whole source bundle, not the flow's own surface state.
                let cold_source = &ResolvedLiquidSurfaceState::try_from_material(
                    cold_law.reference_state().material().clone(),
                    LiquidViscosityBasis::Dynamic,
                    LiquidCapillarity::Neglected,
                )
                .unwrap();
                let before_bits = old.flow().lattice().grid().state_bits();
                let first = old
                    .advance_uniform_slab_volume(
                        face,
                        400_000.0,
                        hot_source,
                        hot_law,
                        controls(),
                        gate,
                        cx,
                    )
                    .unwrap();
                let second = first
                    .liquid
                    .advance_uniform_slab_volume(
                        face,
                        -400_000.0,
                        cold_source,
                        cold_law,
                        controls(),
                        gate,
                        cx,
                    )
                    .unwrap();
                let third = second
                    .liquid
                    .advance_uniform_slab_volume(
                        face,
                        400_000.0,
                        hot_source,
                        hot_law,
                        controls(),
                        gate,
                        cx,
                    )
                    .unwrap();
                let axis = face as usize / 2;
                near(
                    second.report.initial_height_m,
                    first.report.final_height_m,
                    2e-13,
                );
                near(
                    second.report.final_height_m,
                    first.report.initial_height_m,
                    2e-13,
                );
                near(
                    third.report.final_height_m,
                    first.report.final_height_m,
                    2e-13,
                );
                near(
                    second.report.mass_after_kg,
                    first.report.mass_before_kg,
                    2e-13,
                );
                near(
                    third.report.mass_after_kg,
                    first.report.mass_before_kg,
                    2e-13,
                );
                near(
                    second.report.native_kinetic_energy_before_j,
                    first.report.native_kinetic_energy_j,
                    2e-13,
                );
                near(
                    third.report.native_kinetic_energy_before_j,
                    second.report.native_kinetic_energy_j,
                    2e-13,
                );
                for (incoming, outgoing) in [
                    (&first.report, &second.report),
                    (&second.report, &third.report),
                ] {
                    for a in 0..3 {
                        near(
                            outgoing.momentum_before_kg_m_s[a],
                            incoming.momentum_after_kg_m_s[a],
                            2e-12,
                        );
                        assert!(
                            (outgoing.momentum_before_kg_m_s[a]
                                + outgoing.pressure_impulse_kg_m_s[a]
                                - outgoing.momentum_after_kg_m_s[a])
                                .abs()
                                < 1e-16
                        );
                    }
                    assert!(outgoing.mechanical_energy_defect_j.abs() < 1e-18);
                    assert!(outgoing.native_to_face_kinetic_change_j.is_finite());
                }
                assert!(
                    first.report.momentum_after_kg_m_s[axis]
                        * second.report.momentum_after_kg_m_s[axis]
                        < 0.0
                );
                near(
                    second.liquid.total_enthalpy_j().unwrap(),
                    old.total_enthalpy_j().unwrap(),
                    2e-12,
                );
                near(third.liquid.physical_time_s(), 0.253, 1e-14);
                assert_eq!(third.liquid.accepted_steps(), 3);
                assert_eq!(old.flow().lattice().grid().state_bits(), before_bits);
                assert_eq!(first.liquid.accepted_steps(), 1);
            },
        );
    }
}

#[test]
fn g3_zero_source_projects_retained_motion_without_new_mass_heat_or_expansion() {
    with_case(
        Case {
            terminal_fill: 0.55,
            ..Case::default()
        },
        |old, source, law, gate, cx| {
            let first = old
                .advance_uniform_slab_volume(
                    Face3::XMin,
                    400_000.0,
                    source,
                    law,
                    controls(),
                    gate,
                    cx,
                )
                .unwrap();
            let second = first
                .liquid
                .advance_uniform_slab_volume(Face3::XMin, 0.0, source, law, controls(), gate, cx)
                .unwrap();
            near(
                second.report.final_height_m,
                first.report.final_height_m,
                2e-13,
            );
            near(
                second.report.mass_after_kg,
                first.report.mass_after_kg,
                2e-13,
            );
            near(
                second.liquid.total_enthalpy_j().unwrap(),
                first.liquid.total_enthalpy_j().unwrap(),
                2e-13,
            );
            assert_eq!(second.report.requested_enthalpy_increment_j, 0.0);
            assert_eq!(second.report.affine_kinetic_energy_j, 0.0);
            assert_eq!(second.report.native_kinetic_energy_j, 0.0);
            assert_eq!(second.report.required_mechanical_work_j, 0.0);
            assert!(second.projection.projection.report.projection_loss_j > 0.0);
            assert!(second.report.native_kinetic_energy_before_j > 0.0);
            assert!(second.report.mechanical_energy_defect_j.abs() < 1e-18);
            assert!(
                (second.report.pressure_impulse_kg_m_s[0]
                    + second.report.momentum_before_kg_m_s[0])
                    .abs()
                    < 1e-16
            );
            let no_mapping = FreeSurfaceSlabVolumeControls {
                maximum_kinetic_mapping_change_j: 0.0,
                ..controls()
            };
            assert!(
                first
                    .liquid
                    .advance_uniform_slab_volume(
                        Face3::XMin,
                        0.0,
                        source,
                        law,
                        no_mapping,
                        gate,
                        cx,
                    )
                    .is_err()
            );
            let opposite = first.liquid.advance_uniform_slab_volume(
                Face3::XMax,
                0.0,
                source,
                law,
                controls(),
                gate,
                cx,
            );
            assert!(opposite.is_err());
            let mut evolved = first.liquid.clone();
            evolved.step(gate).unwrap();
            assert!(
                evolved
                    .advance_uniform_slab_volume(
                        Face3::XMin,
                        0.0,
                        source,
                        law,
                        controls(),
                        gate,
                        cx,
                    )
                    .is_err()
            );
            assert_eq!(first.liquid.accepted_steps(), 1);
        },
    );
}
