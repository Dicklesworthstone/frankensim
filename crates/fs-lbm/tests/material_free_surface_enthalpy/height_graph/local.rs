//! G1/G4: detached-drop stress consumes the public sourced material/VOF owner.
//! The manufactured sphere is geometry evidence, not an experimental dataset.

use super::*;
use fs_lbm::material::free_surface::enthalpy::thermocapillary::LocalThermocapillarity;

#[cfg(feature = "enthalpy-coupling")]
#[path = "local/surface_heat.rs"]
mod surface_heat;

const SIZE: usize = 56;
const CENTER: usize = SIZE / 2;
const RADIUS: f64 = 16.0;
const TANGENT_POLE: [usize; 3] = [CENTER, CENTER, CENTER + 16];
const NORMAL_POLE: [usize; 3] = [CENTER + 16, CENTER, CENTER];

fn local() -> LocalThermocapillarity {
    LocalThermocapillarity {
        minimum_surface_fill: 1e-10,
        maximum_gradient_condition: 1e4,
    }
}

fn drop_domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [SIZE; 3],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: SIZE.pow(3),
    }
}

// Independent quadrature of the exact sphere's vertical liquid interval. This
// fixture does not call the native normal, plane/cube-cut or traction routines.
fn sphere_fill(coordinates: [usize; 3]) -> f64 {
    let p = coordinates.map(|v| v as f64 - CENTER as f64);
    let square = RADIUS * RADIUS;
    if p.iter()
        .map(|x| (x.abs() - 0.5).max(0.0).powi(2))
        .sum::<f64>()
        >= square
    {
        return 0.0;
    }
    if p.iter().map(|x| (x.abs() + 0.5).powi(2)).sum::<f64>() <= square {
        return 1.0;
    }
    const NODES: [f64; 4] = [
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526,
    ];
    const WEIGHTS: [f64; 4] = [
        0.3478548451374538,
        0.6521451548625461,
        0.6521451548625461,
        0.3478548451374538,
    ];
    let mut volume = 0.0;
    // Subdivide before Gaussian integration to resolve a clipped surface near
    // a cell face; unsplit nodes can spuriously give an exactly half-full pole.
    for i in 0..32 {
        for j in 0..32 {
            let (a, b) = (i % 4, j % 4);
            let u = p[0] - 0.5 + (i / 4) as f64 / 8.0 + (0.5 + 0.5 * NODES[a]) / 8.0;
            let v = p[1] - 0.5 + (j / 4) as f64 / 8.0 + (0.5 + 0.5 * NODES[b]) / 8.0;
            let radial = square - u * u - v * v;
            if radial > 0.0 {
                let height = radial.sqrt();
                let interval = ((p[2] + 0.5).min(height) - (p[2] - 0.5).max(-height)).max(0.0);
                volume += WEIGHTS[a] * WEIGHTS[b] * interval / 256.0;
            }
        }
    }
    volume.min(1.0)
}

fn drop_cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut cells = Vec::new();
    for z in CENTER - 16..=CENTER + 16 {
        for y in CENTER - 16..=CENTER + 16 {
            for x in CENTER - 16..=CENTER + 16 {
                let fill_fraction = sphere_fill([x, y, z]);
                if fill_fraction > 0.0 {
                    cells.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    cells
}

fn with_drop<R>(
    slope: f64,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    with_moving_drop(slope, [0.0; 3], [1.0, 0.0, 0.0], run)
}

fn with_moving_drop<R>(
    slope: f64,
    velocity_m_s: [f64; 3],
    temperature_gradient_k: [f64; 3],
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = source_card(slope);
    let source = source(&card);
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let caloric = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut cells = drop_cells();
        for cell in &mut cells {
            cell.velocity_m_s = velocity_m_s;
        }
        let mut thermal = plan_free_surface_from_material(&source, DX, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_local_thermocapillarity(
                drop_domain(),
                &cells,
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls {
                    max_cells: 65536,
                    max_mass_routes: 1_000_000,
                    ..FreeSurfaceEnthalpyControls::default()
                },
                local(),
                None,
                None,
                |p| {
                    200000.0
                        + 2000.0
                            * (0..3)
                                .map(|d| temperature_gradient_k[d] * (p[d] as f64 - CENTER as f64))
                                .sum::<f64>()
                },
                cx,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

#[test]
fn g1_local_drop_uses_sourced_tangent_stress_and_surface_centroid_pressure() {
    let slope = 1e-6;
    with_drop(slope, |thermal, gate| {
        assert_eq!(thermal.local_thermocapillary_surface(), Some(local()));
        assert!(thermal.thermocapillary_surface().is_none());
        assert!(thermal.height_graph_thermocapillary_surface().is_none());

        // Symmetry makes the north-pole normal exactly z. The ambient physical
        // gradient is slope/dx; a unit projected patch loads its actual liquid
        // volume. This checks SI conversion and removes dependence on a guessed
        // spherical area or on the production surface-gradient implementation.
        let [x, y, z] = TANGENT_POLE;
        let force = thermal.flow().lattice().force_density_at(x, y, z).unwrap();
        let expected = slope * DT * DT / (1000.0 * DX.powi(3) * sphere_fill(TANGENT_POLE));
        near(force[0], expected, 2e-10);
        assert!(force[1].abs() < 1e-16 && force[2].abs() < 1e-16);
        assert!(
            thermal
                .flow()
                .sample_at_si(x, y, z)
                .unwrap()
                .unwrap()
                .velocity_m_s
                .iter()
                .all(|u| u.abs() < 1e-15),
            "binding a stress must preserve the specified initial velocity"
        );

        // At the east pole the gradient is purely normal: it changes the
        // extrapolated surface pressure coefficient but supplies no Marangoni
        // traction. Bulk source receipts stay at their actual caloric T=416 K.
        let [x, y, z] = NORMAL_POLE;
        assert!(
            thermal
                .flow()
                .lattice()
                .force_density_at(x, y, z)
                .unwrap()
                .iter()
                .all(|f| f.abs() < 1e-16)
        );
        let bulk = thermal
            .material_at(NORMAL_POLE)
            .unwrap()
            .material()
            .property(SURFACE_TENSION_PROPERTY)
            .unwrap();
        near(bulk.value_si(), 0.001 + 16.0 * slope, 2e-14);
        let pin = bulk.answer().receipt.selected;
        let surface = thermal
            .flow()
            .lattice()
            .surface_tension_at(x, y, z)
            .unwrap();
        let expected = (bulk.value_si() + slope * (sphere_fill(NORMAL_POLE) - 0.5)) * DT * DT
            / (1000.0 * DX.powi(3));
        near(surface, expected, 2e-12);
        assert!((surface - bulk.value_si()).abs() > 1e-9);

        let mass = thermal.flow().ledger_mass_kg().unwrap();
        let energy = thermal.total_enthalpy_j().unwrap();
        let report = thermal.step(gate).unwrap();
        near(report.mass_after_kg, mass, 3e-12);
        near(report.enthalpy_after_j, energy, 3e-11);
        assert_eq!(report.external_heat_j, 0.0);
        assert_eq!(report.viscous_heat_j, 0.0);
        assert!(report.balance_defect_j.abs() < 1e-7);
        assert_eq!(
            thermal
                .material_at(NORMAL_POLE)
                .unwrap()
                .material()
                .property(SURFACE_TENSION_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected,
            pin
        );
        let [x, y, z] = TANGENT_POLE;
        assert!(
            thermal
                .flow()
                .sample_at_si(x, y, z)
                .unwrap()
                .unwrap()
                .velocity_m_s[0]
                > 0.0,
            "the original flow owner must move toward higher surface tension"
        );
    });
}

#[test]
fn g4_local_drop_source_refusal_and_cancellation_leave_accepted_fields_unchanged() {
    with_drop(1e-6, |thermal, _| {
        let before = thermal.flow().lattice().grid().state_bits();
        let mass = thermal.flow().ledger_mass_kg().unwrap().to_bits();
        let energy = thermal.total_enthalpy_j().unwrap().to_bits();
        let state = thermal.material_at(TANGENT_POLE).unwrap().clone();
        let enthalpy = thermal
            .specific_enthalpy_at(TANGENT_POLE)
            .unwrap()
            .to_bits();
        let [x, y, z] = TANGENT_POLE;
        let force = thermal.flow().lattice().force_density_at(x, y, z);
        let sigma = thermal.flow().lattice().surface_tension_at(x, y, z);
        let gate = CancelGate::new_clock_free();
        // T=550 K lies inside the caloric chart but outside every source claim.
        // Failure must therefore come from material support, not a bad heat input.
        let heat = sphere_fill(TANGENT_POLE) * 1000.0 * DX.powi(3) * 300000.0;
        let result = thermal.step_with_heat(
            &[FreeSurfaceCellHeat {
                coordinates: TANGENT_POLE,
                energy_j: heat,
            }],
            &gate,
        );
        assert!(
            matches!(result, Err(FreeSurfaceEnthalpyError::Material(_))),
            "{result:?}"
        );
        for cancel in [false, true] {
            if cancel {
                gate.request();
                assert!(thermal.step(&gate).is_err());
            }
            assert_eq!(before, thermal.flow().lattice().grid().state_bits());
            assert_eq!(mass, thermal.flow().ledger_mass_kg().unwrap().to_bits());
            assert_eq!(energy, thermal.total_enthalpy_j().unwrap().to_bits());
            assert_eq!(&state, thermal.material_at(TANGENT_POLE).unwrap());
            assert_eq!(
                enthalpy,
                thermal
                    .specific_enthalpy_at(TANGENT_POLE)
                    .unwrap()
                    .to_bits()
            );
            assert_eq!(force, thermal.flow().lattice().force_density_at(x, y, z));
            assert_eq!(sigma, thermal.flow().lattice().surface_tension_at(x, y, z));
            assert_eq!(thermal.local_thermocapillary_surface(), Some(local()));
        }
    });
}

#[test]
fn g4_local_constructor_refuses_invalid_resolution_controls_without_changing_sources() {
    let card = source_card(1e-6);
    let source = source(&card);
    let before = source.material().clone();
    let law = LiquidTemperatureLaw::try_from_material(
        &card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let caloric = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, DX, DT, 0.1).unwrap();
        let initial = graph_cells(false);
        // Invalid controls are diagnosed before local reconstruction, so this
        // inexpensive fixture can exercise every admission boundary.
        for (minimum_surface_fill, maximum_gradient_condition) in [
            (0.0, 1e4),
            (1.0, 1e4),
            (f64::NAN, 1e4),
            (1e-10, 0.5),
            (1e-10, f64::INFINITY),
        ] {
            let result = scaling.free_surface3_from_cells_with_local_thermocapillarity(
                domain(),
                &initial,
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                LocalThermocapillarity {
                    minimum_surface_fill,
                    maximum_gradient_condition,
                },
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            );
            let error = match result {
                Ok(_) => panic!("invalid local reconstruction controls were accepted"),
                Err(error) => error,
            };
            assert!(
                error
                    .to_string()
                    .contains("bounded reconstruction conditioning"),
                "{error}"
            );
            assert_eq!(source.material(), &before);
        }
    });
}
