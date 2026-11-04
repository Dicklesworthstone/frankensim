//! G1/G3: physical outlet flow, coordinate orientation and complete storage.
use fs_math::c64::C64;
use fs_phs::acoustic_outlet::{
    AcousticCellWall, AcousticSectionMedium, AcousticWallRelaxation, acoustic_cell_geometry,
    acoustic_chain_with_outlet_flow, acoustic_chain_with_outlet_flow_in_exterior,
    acoustic_chain_with_regional_outlet_flow,
    acoustic_chain_with_regional_outlet_flow_and_cell_walls,
};
use fs_phs::harmonic::{HarmonicBudget, HarmonicModel};
use fs_phs::{
    AcousticSection, AcousticTap, ChainOutlet, MouthFlange, ViscothermalPin, WallPin,
    acoustic_chain_outlet_wall,
};

fn section(cells: usize, cone: bool) -> AcousticSection {
    AcousticSection {
        length: 0.3,
        radius: 0.01,
        outlet_radius: if cone { 0.02 } else { 0.01 },
        cells,
    }
}

#[test]
fn long_open_chimney_has_an_open_low_frequency_limit_and_conserves_total_volume() {
    let density = 1.2;
    let speed = 343.0;
    let radius = 0.004;
    for cells in [3, 4] {
        let sections = [AcousticSection {
            length: 0.34,
            radius: 0.012,
            outlet_radius: 0.012,
            cells,
        }];
        let tap = AcousticTap {
            station: 0.6,
            neck_length: 0.05,
            neck_radius: radius,
            open_fraction: 1.0,
        };
        let (system, bell) = acoustic_chain_with_outlet_flow(
            &sections,
            density,
            speed,
            ChainOutlet::Rigid,
            1,
            &[tap],
            None,
            None,
        )
        .unwrap();
        let base = 2 * cells;
        let outlet_flux = base + 4;
        let tap_cell = (tap.station * cells as f64).floor() as usize;
        let sign = if tap_cell % 2 == 0 { 1.0 } else { -1.0 };
        assert_eq!(system.state_dim(), base + 5);
        let model = HarmonicModel::from_ode(
            &system,
            HarmonicBudget {
                max_dimension: 64,
                max_dense_bytes: 2_000_000,
                max_dense_work: 2_000_000,
                relative_residual_tolerance: 1e-10,
            },
            || false,
        )
        .unwrap();
        let omega_pin = core::f64::consts::PI * speed / (2.0 * sections[0].length);
        // The original compact radiation owner retains its fixed pin. Its
        // finite R is the only real load in this otherwise lossless fixture.
        let resistance = density * omega_pin * omega_pin / (2.0 * core::f64::consts::PI * speed);
        let low = model
            .factor(0.001, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap();
        let outward = low.effort[outlet_flux].scale(sign);
        assert!((outward - C64::ONE).abs() < 1e-6);
        assert!((low.output[0].re / resistance - 1.0).abs() < 1e-6);
        assert!(
            low.output[0].im.abs() < 0.01 * resistance,
            "an open chimney cannot retain a divergent sealed-cap compliance"
        );

        let omega = 700.0;
        let wave = model
            .factor(omega, || false)
            .unwrap()
            .response(&[C64::ONE], || false)
            .unwrap();
        let neck = wave.effort[base + 1].scale(sign);
        let outward = wave.effort[outlet_flux].scale(sign);
        let compression = (wave.state[base] + wave.state[base + 2]).scale(sign);
        assert!(
            (neck - outward - C64::new(0.0, -omega) * compression).abs()
                < 1e-9 * (neck.abs() + outward.abs())
        );
        assert!(
            (neck - outward).abs() > 1e-5 * neck.abs(),
            "distributed chimney compression makes inlet flow different from mouth flow"
        );

        let mut initial = vec![0.0; system.state_dim()];
        for cell in 0..cells {
            let orientation = if cell % 2 == 0 { 1.0 } else { -1.0 };
            initial[2 * cell] = orientation * (cell + 1) as f64 * 1e-8;
        }
        initial[base] = sign * 2e-10;
        initial[base + 2] = sign * 3e-10;
        initial[outlet_flux] = sign * 1e-4;
        let dt = 1e-6;
        let input = 3e-4;
        let step = fs_phs::step(&system, &initial, &[input], dt).unwrap();
        let midpoint: Vec<_> = initial
            .iter()
            .zip(&step.x)
            .map(|(a, b)| f64::midpoint(*a, *b))
            .collect();
        let outward = sign * system.effort(&midpoint)[outlet_flux];
        let mut compression_change =
            sign * (step.x[base] - initial[base] + step.x[base + 2] - initial[base + 2]);
        for cell in 0..cells {
            let orientation = if cell % 2 == 0 { 1.0 } else { -1.0 };
            compression_change += orientation * (step.x[2 * cell] - initial[2 * cell]);
        }
        assert_eq!(bell.volume_flow_m3_s(&midpoint).unwrap(), 0.0);
        let boundary = dt * (input - outward);
        assert!((compression_change - boundary).abs() < 1e-10 * boundary.abs().max(1e-12));
        let radiation_work = dt * resistance * outward * outward;
        assert!((step.dissipated - radiation_work).abs() < 1e-10 * radiation_work.max(1e-20));
    }
}

#[test]
fn mixed_chimney_order_preserves_response_with_all_original_memory_and_wall_tails() {
    let short = AcousticTap {
        station: 0.3,
        neck_length: 0.003,
        neck_radius: 0.003,
        open_fraction: 0.5,
    };
    let long = AcousticTap {
        station: 0.7,
        neck_length: 0.05,
        neck_radius: 0.004,
        open_fraction: 1.0,
    };
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 2,
    };
    let wall = WallPin {
        surface_density: 0.1,
        stiffness_per_area: 1e7,
        resistance: 10.0,
    };
    for cone in [false, true] {
        // The single cone exercises the original spherical caller as well as
        // the cylindrical frustum caller; both share the corrected chimney.
        let sections = [section(4, cone)];
        let response = |taps: &[AcousticTap]| {
            let system = acoustic_chain_outlet_wall(
                &sections,
                1.2,
                343.0,
                ChainOutlet::Radiating(MouthFlange::Flanged),
                1,
                taps,
                Some(&pin),
                Some(&wall),
            )
            .unwrap();
            let model = HarmonicModel::from_ode(
                &system,
                HarmonicBudget {
                    max_dimension: 128,
                    max_dense_bytes: 8_000_000,
                    max_dense_work: 20_000_000,
                    relative_residual_tolerance: 1e-10,
                },
                || false,
            )
            .unwrap();
            model
                .factor(900.0, || false)
                .unwrap()
                .response(&[C64::ONE], || false)
                .unwrap()
        };
        let forward = response(&[short, long]);
        let reversed = response(&[long, short]);
        assert!((forward.output[0] - reversed.output[0]).abs() < 1e-9 * forward.output[0].abs());
        assert!(forward.dissipated_power_w > 0.0);
        assert!(
            (forward.dissipated_power_w - reversed.dissipated_power_w).abs()
                < 1e-9 * forward.dissipated_power_w
        );
    }
}

#[test]
fn source_wall_geometry_uses_surface_slant_without_changing_axial_gas_storage() {
    let pin = WallPin {
        surface_density: 0.08,
        stiffness_per_area: 2e6,
        resistance: 12.0,
    };
    let medium = AcousticSectionMedium {
        density: 1.2,
        sound_speed: 340.0,
        viscothermal: None,
    };
    for cone in [false, true] {
        let sections = [AcousticSection {
            length: 0.08,
            radius: 0.01,
            outlet_radius: if cone { 0.03 } else { 0.01 },
            cells: 3,
        }];
        let geometry = acoustic_cell_geometry(&sections).unwrap();
        let walls = vec![
            AcousticCellWall {
                wall: pin,
                relaxation: vec![]
            };
            3
        ];
        let (system, flow) = acoustic_chain_with_regional_outlet_flow_and_cell_walls(
            &sections,
            &[medium],
            1.2,
            340.0,
            ChainOutlet::PressureRelease,
            1,
            &walls,
        )
        .unwrap();
        let (legacy, legacy_flow) = acoustic_chain_with_regional_outlet_flow(
            &sections,
            &[medium],
            1.2,
            340.0,
            ChainOutlet::PressureRelease,
            1,
            Some(&pin),
        )
        .unwrap();
        assert_eq!(system.state_dim(), 12);
        let diagonal = system.effort(&[1.0; 12]);
        let old = legacy.effort(&[1.0; 12]);
        let dx = sections[0].length / 3.0;
        let dr = (sections[0].outlet_radius - sections[0].radius) / 3.0;
        for (cell, geometry) in geometry.iter().enumerate() {
            let ra = sections[0].radius + cell as f64 * dr;
            let rb = ra + dr;
            let area = core::f64::consts::PI * (ra + rb) * dx.hypot(dr);
            assert_eq!(geometry.section, 0);
            assert_eq!(geometry.axial_length.to_bits(), dx.to_bits());
            assert!((geometry.wall_area_m2() / area - 1.0).abs() < 1e-14);
            assert!((diagonal[6 + 2 * cell] * area / pin.stiffness_per_area - 1.0).abs() < 1e-13);
            assert!((diagonal[7 + 2 * cell] * pin.surface_density / area - 1.0).abs() < 1e-13);
            // Surface slant cannot alter the gas volume or axial inertia.
            assert_eq!(diagonal[2 * cell].to_bits(), old[2 * cell].to_bits());
            assert_eq!(
                diagonal[2 * cell + 1].to_bits(),
                old[2 * cell + 1].to_bits()
            );
        }
        if !cone {
            let (j, r, g) = system.structure();
            let (lj, lr, lg) = legacy.structure();
            let bits = |values: &[f64]| values.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(j), bits(lj));
            assert_eq!(bits(r), bits(lr));
            assert_eq!(bits(g), bits(lg));
            assert_eq!(bits(&diagonal), bits(&old));
            assert_eq!(bits(flow.coefficients()), bits(legacy_flow.coefficients()));
        } else {
            assert!(diagonal[6] < old[6]);
        }
    }
}

#[test]
fn cell_wall_maxwell_pressure_work_and_flow_match_the_complete_native_state() {
    for cells in [3, 4] {
        let sections = [section(cells, true)];
        let geometry = acoustic_cell_geometry(&sections).unwrap();
        let medium = AcousticSectionMedium {
            density: 1.2,
            sound_speed: 340.0,
            viscothermal: None,
        };
        let walls: Vec<_> = (0..cells)
            .map(|cell| AcousticCellWall {
                wall: WallPin {
                    surface_density: 0.05,
                    stiffness_per_area: 1e6,
                    resistance: 10.0,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: (cell + 1) as f64 * 2e5,
                    relaxation_time_s: 0.002,
                }],
            })
            .collect();
        let build = |walls: &[AcousticCellWall]| {
            acoustic_chain_with_regional_outlet_flow_and_cell_walls(
                &sections,
                &[medium],
                1.2,
                340.0,
                ChainOutlet::PressureRelease,
                1,
                walls,
            )
        };
        let (system, flow) = build(&walls).unwrap();
        assert_eq!(system.state_dim(), 5 * cells);
        assert_eq!(flow.coefficients().len(), system.state_dim());
        let model = HarmonicModel::from_ode(
            &system,
            HarmonicBudget {
                max_dimension: 64,
                max_dense_bytes: 2_000_000,
                max_dense_work: 2_000_000,
                relative_residual_tolerance: 1e-10,
            },
            || false,
        )
        .unwrap();
        for omega in [200.0, 1500.0, 6000.0] {
            let response = model
                .factor(omega, || false)
                .unwrap()
                .response(&[C64::ONE], || false)
                .unwrap();
            let s = C64::new(0.0, -omega);
            for (cell, wall) in walls.iter().enumerate() {
                let area = geometry[cell].wall_area_m2();
                let arm = wall.relaxation[0];
                // Independent local mechanical impedance in pressure/volume-flow
                // units under exp(-i omega t); no fitted wall admittance.
                let impedance = C64::new(
                    wall.wall.resistance / area,
                    -omega * wall.wall.surface_density / area,
                ) + C64::new(wall.wall.stiffness_per_area / area, 0.0) / s
                    + C64::new(arm.stiffness_per_area / area, 0.0)
                        / (s + C64::new(1.0 / arm.relaxation_time_s, 0.0));
                let mut basis = vec![0.0; system.state_dim()];
                basis[2 * cell] = 1.0;
                let pressure = response.state[2 * cell].scale(system.effort(&basis)[2 * cell]);
                let volume_flow = s * response.state[2 * cells + 2 * cell];
                assert!((pressure - impedance * volume_flow).abs() < 1e-9 * pressure.abs());
            }
        }
        let mut initial = vec![0.0; system.state_dim()];
        for cell in 0..cells {
            initial[2 * cell] = 1e-8 * (cell + 1) as f64;
            initial[2 * cells + 2 * cell] = 1e-9;
            initial[4 * cells + cell] = 1e-4;
        }
        let dt = 1e-6;
        let step = fs_phs::step(&system, &initial, &[3e-4], dt).unwrap();
        let midpoint: Vec<_> = initial
            .iter()
            .zip(&step.x)
            .map(|(a, b)| f64::midpoint(*a, *b))
            .collect();
        let mut expected_loss = 0.0;
        let mut volume_change = 0.0;
        for (cell, wall) in walls.iter().enumerate() {
            let area = geometry[cell].wall_area_m2();
            let q_wall = 2 * cells + 2 * cell;
            let velocity = midpoint[q_wall + 1] * area / wall.wall.surface_density;
            let arm = wall.relaxation[0];
            let elastic_memory = (arm.stiffness_per_area / area).sqrt() * midpoint[q_wall]
                - midpoint[4 * cells + cell];
            expected_loss += dt
                * (wall.wall.resistance / area * velocity.powi(2)
                    + elastic_memory.powi(2) / arm.relaxation_time_s);
            let sign = if cell % 2 == 0 { 1.0 } else { -1.0 };
            volume_change +=
                sign * (step.x[2 * cell] - initial[2 * cell] + step.x[q_wall] - initial[q_wall]);
        }
        assert!((step.dissipated - expected_loss).abs() < 1e-10 * expected_loss);
        let boundary = dt * (3e-4 - flow.volume_flow_m3_s(&midpoint).unwrap());
        assert!((volume_change - boundary).abs() < 1e-10 * boundary.abs().max(1e-12));
        let mut invalid = walls.clone();
        invalid[0].relaxation[0].relaxation_time_s = 0.0;
        assert!(build(&invalid).is_err());
        assert!(build(&walls[..cells - 1]).is_err());
    }
}

fn geometries(cells: usize) -> [(&'static str, Vec<AcousticSection>); 3] {
    [
        ("cylinder", vec![section(cells, false)]),
        ("single cone", vec![section(cells, true)]),
        (
            "joined frusta",
            vec![
                AcousticSection {
                    length: 0.1,
                    radius: 0.01,
                    outlet_radius: 0.015,
                    cells: 1,
                },
                AcousticSection {
                    length: 0.2,
                    radius: 0.015,
                    outlet_radius: 0.02,
                    cells: cells - 1,
                },
            ],
        ),
    ]
}

#[test]
fn outward_flow_closes_independent_compression_volume_balance_for_odd_and_even_chains() {
    for cells in [2, 3, 4, 5] {
        for (shape, sections) in geometries(cells) {
            for outlet in [
                ChainOutlet::PressureRelease,
                ChainOutlet::Radiating(MouthFlange::Flanged),
            ] {
                let (system, flow) = acoustic_chain_with_outlet_flow(
                    &sections,
                    1.2,
                    340.0,
                    outlet,
                    1,
                    &[],
                    None,
                    None,
                )
                .unwrap();
                let mut initial = vec![0.0; system.state_dim()];
                for i in 0..cells {
                    initial[2 * i] = (i as f64 + 1.0) * 1e-8;
                    initial[2 * i + 1] = (2 * i + 1) as f64 * 1e-4;
                }
                let input = 2e-4;
                let dt = 1e-5;
                let step = fs_phs::step(&system, &initial, &[input], dt).unwrap();
                let midpoint: Vec<_> = initial
                    .iter()
                    .zip(&step.x)
                    .map(|(a, b)| f64::midpoint(*a, *b))
                    .collect();
                let volume_change: f64 = (0..cells)
                    .map(|i| {
                        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                        sign * (step.x[2 * i] - initial[2 * i])
                    })
                    .sum();
                let boundary_change = dt * (input - flow.volume_flow_m3_s(&midpoint).unwrap());
                assert!(
                    (volume_change - boundary_change).abs()
                        < 1e-11 * boundary_change.abs().max(1e-12),
                    "{shape}, cells={cells}, {outlet:?}: stored={volume_change:e}, boundary={boundary_change:e}"
                );
            }
        }
    }
}

#[test]
fn cylinder_and_frustum_preserve_physical_low_frequency_through_flow() {
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 2,
    };
    for cells in [3, 4] {
        for (shape, sections) in geometries(cells) {
            for loss in [None, Some(&pin)] {
                let (system, flow) = acoustic_chain_with_outlet_flow(
                    &sections,
                    1.2,
                    340.0,
                    ChainOutlet::Radiating(MouthFlange::Flanged),
                    1,
                    &[],
                    loss,
                    None,
                )
                .unwrap();
                let model = HarmonicModel::from_ode(
                    &system,
                    HarmonicBudget {
                        max_dimension: 256,
                        max_dense_bytes: 32_000_000,
                        max_dense_work: 100_000_000,
                        relative_residual_tolerance: 1e-10,
                    },
                    || false,
                )
                .unwrap();
                let response = model
                    .factor(0.01, || false)
                    .unwrap()
                    .response(&[C64::ONE], || false)
                    .unwrap();
                let outward = flow
                    .coefficients()
                    .iter()
                    .zip(&response.state)
                    .fold(C64::ZERO, |sum, (a, x)| sum + x.scale(*a));
                // Independently, every nonleaking duct transports its inlet
                // flow to the outlet as omega -> 0, regardless of taper or
                // the arbitrary alternating orientation of its native states.
                assert!(
                    (outward - C64::ONE).abs() < 2e-5,
                    "{shape}, cells={cells}, loss={}: {outward:?}",
                    loss.is_some()
                );
                let real: Vec<_> = response.state.iter().map(|x| x.re).collect();
                assert!((flow.volume_flow_m3_s(&real).unwrap() - outward.re).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn observed_terminal_flow_reproduces_independent_physical_radiation_work() {
    let density = 1.2;
    let speed = 340.0;
    let dt = 1e-5;
    for cells in [3, 4] {
        for (shape, sections) in geometries(cells) {
            for (exterior_density, exterior_speed) in [(density, speed), (1.9, 280.0)] {
                let (system, flow) = acoustic_chain_with_outlet_flow_in_exterior(
                    &sections,
                    density,
                    speed,
                    exterior_density,
                    exterior_speed,
                    ChainOutlet::Radiating(MouthFlange::Flanged),
                    1,
                    &[],
                    None,
                    None,
                )
                .unwrap();
                let radius = sections.last().unwrap().outlet_radius;
                let length: f64 = sections.iter().map(|section| section.length).sum();
                // Independently, the compact flanged load at the original pin
                // frequency dissipates R * U_out^2. The accepted step's original
                // material resistance account must equal this physical work.
                let omega = core::f64::consts::PI * speed / (2.0 * length);
                let area = core::f64::consts::PI * radius * radius;
                let resistance = exterior_density * exterior_speed / area
                    * (omega * radius / exterior_speed).powi(2)
                    / 2.0;
                // Every cell keeps its actual interior frustum volume/inertance.
                // Only the terminal flux acquires the original compact mouth mass
                // M_rad = rho_exterior * 0.8216 * radius / area.
                let effort = system.effort(&vec![1.0; system.state_dim()]);
                let mut cell = 0;
                for section in &sections {
                    let dx = section.length / section.cells as f64;
                    for i in 0..section.cells {
                        let ra = section.radius
                            + (section.outlet_radius - section.radius)
                                * (i as f64 / section.cells as f64);
                        let rb = section.radius
                            + (section.outlet_radius - section.radius)
                                * ((i + 1) as f64 / section.cells as f64);
                        let volume =
                            core::f64::consts::PI * (ra * ra + ra * rb + rb * rb) * dx / 3.0;
                        let mut inertia = density * dx / (core::f64::consts::PI * ra * rb);
                        if cell + 1 == cells {
                            inertia += exterior_density * 0.8216 * radius / area;
                        }
                        assert!(
                            (effort[2 * cell] * volume / (density * speed * speed) - 1.0).abs()
                                < 1e-13
                        );
                        assert!((effort[2 * cell + 1] * inertia - 1.0).abs() < 1e-13);
                        cell += 1;
                    }
                }
                let mut initial = vec![0.0; system.state_dim()];
                for i in 0..cells {
                    initial[2 * i] = (i as f64 + 1.0) * 1e-8;
                    initial[2 * i + 1] = (2 * i + 1) as f64 * 1e-4;
                }
                let step = fs_phs::step(&system, &initial, &[2e-4], dt).unwrap();
                let midpoint: Vec<_> = initial
                    .iter()
                    .zip(&step.x)
                    .map(|(a, b)| f64::midpoint(*a, *b))
                    .collect();
                let outward = flow.volume_flow_m3_s(&midpoint).unwrap();
                let expected = dt * resistance * outward * outward;
                let stored_volume_change: f64 = (0..cells)
                    .map(|i| {
                        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                        sign * (step.x[2 * i] - initial[2 * i])
                    })
                    .sum();
                let boundary_change = dt * (2e-4 - outward);
                assert!(
                    (stored_volume_change - boundary_change).abs()
                        < 1e-11 * boundary_change.abs().max(1e-12)
                );
                assert!(expected > 0.0);
                assert!(
                    (step.dissipated - expected).abs() < 1e-11 * expected,
                    "{shape}, cells={cells}: native={}, physical={expected}",
                    step.dissipated
                );
            }
        }
    }
}

#[test]
fn explicitly_identical_exterior_preserves_native_storage_structure_and_flow_bits() {
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 2,
    };
    let bits = |values: &[f64]| {
        values
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    for (_, sections) in geometries(3) {
        for loss in [None, Some(&pin)] {
            let (original, original_flow) = acoustic_chain_with_outlet_flow(
                &sections,
                1.2,
                340.0,
                ChainOutlet::Radiating(MouthFlange::Flanged),
                1,
                &[],
                loss,
                None,
            )
            .unwrap();
            let (explicit, explicit_flow) = acoustic_chain_with_outlet_flow_in_exterior(
                &sections,
                1.2,
                340.0,
                1.2,
                340.0,
                ChainOutlet::Radiating(MouthFlange::Flanged),
                1,
                &[],
                loss,
                None,
            )
            .unwrap();
            let media = vec![
                AcousticSectionMedium {
                    density: 1.2,
                    sound_speed: 340.0,
                    viscothermal: loss.copied(),
                };
                sections.len()
            ];
            let (regional, regional_flow) = acoustic_chain_with_regional_outlet_flow(
                &sections,
                &media,
                1.2,
                340.0,
                ChainOutlet::Radiating(MouthFlange::Flanged),
                1,
                None,
            )
            .unwrap();
            let (j, r, g) = original.structure();
            let (ej, er, eg) = explicit.structure();
            let (rj, rr, rg) = regional.structure();
            assert_eq!(bits(j), bits(ej));
            assert_eq!(bits(r), bits(er));
            assert_eq!(bits(g), bits(eg));
            assert_eq!(bits(j), bits(rj));
            assert_eq!(bits(r), bits(rr));
            assert_eq!(bits(g), bits(rg));
            assert_eq!(
                bits(original_flow.coefficients()),
                bits(explicit_flow.coefficients())
            );
            assert_eq!(
                bits(original_flow.coefficients()),
                bits(regional_flow.coefficients())
            );
            for i in 0..original.state_dim() {
                let mut basis = vec![0.0; original.state_dim()];
                basis[i] = 1.0;
                assert_eq!(
                    bits(&original.effort(&basis)),
                    bits(&explicit.effort(&basis))
                );
                assert_eq!(
                    bits(&original.effort(&basis)),
                    bits(&regional.effort(&basis))
                );
            }
        }
    }
}

#[test]
fn regional_frusta_keep_local_integrals_interface_power_and_exterior_radiation_work() {
    for last_cells in [2, 3] {
        let sections = [
            AcousticSection {
                length: 0.1,
                radius: 0.01,
                outlet_radius: 0.015,
                cells: 2,
            },
            AcousticSection {
                length: 0.2,
                radius: 0.015,
                outlet_radius: 0.02,
                cells: last_cells,
            },
        ];
        let media = [
            AcousticSectionMedium {
                density: 0.8,
                sound_speed: 420.0,
                viscothermal: None,
            },
            AcousticSectionMedium {
                density: 2.1,
                sound_speed: 270.0,
                viscothermal: None,
            },
        ];
        let (system, flow) = acoustic_chain_with_regional_outlet_flow(
            &sections,
            &media,
            1.2,
            340.0,
            ChainOutlet::Radiating(MouthFlange::Flanged),
            1,
            None,
        )
        .unwrap();
        let cells = 2 + last_cells;
        let n = system.state_dim();
        assert_eq!(n, 2 * cells);
        let mut initial = vec![0.0; n];
        let coefficients = system.effort(&vec![1.0; n]);
        let mut cell = 0;
        for (section, medium) in sections.iter().zip(media) {
            let dx = section.length / section.cells as f64;
            for i in 0..section.cells {
                let dr = section.outlet_radius - section.radius;
                let ra = section.radius + dr * (i as f64 / section.cells as f64);
                let rb = section.radius + dr * ((i + 1) as f64 / section.cells as f64);
                let volume = core::f64::consts::PI * (ra * ra + ra * rb + rb * rb) * dx / 3.0;
                let compliance = volume / (medium.density * medium.sound_speed.powi(2));
                let mut inertia = medium.density * dx / (core::f64::consts::PI * ra * rb);
                if cell + 1 == cells {
                    inertia += 1.2 * 0.8216 / (core::f64::consts::PI * rb);
                }
                assert!((coefficients[2 * cell] * compliance - 1.0).abs() < 1e-13);
                assert!((coefficients[2 * cell + 1] * inertia - 1.0).abs() < 1e-13);
                let sign = if cell % 2 == 0 { 1.0 } else { -1.0 };
                initial[2 * cell] = sign * compliance * (100.0 + 7.0 * cell as f64);
                initial[2 * cell + 1] = -sign * inertia * 2e-4;
                cell += 1;
            }
        }
        let effort = system.effort(&initial);
        let (j, _, _) = system.structure();
        // The final series branch of region 0 meets the SAME pressure node
        // as region 1. Its physical momentum changes by p_left - p_right,
        // without any density-ratio multiplier or added interface inertance.
        let interface_flux = 3;
        let momentum_rate: f64 = j[interface_flux * n..(interface_flux + 1) * n]
            .iter()
            .zip(&effort)
            .map(|(a, b)| a * b)
            .sum();
        assert!((momentum_rate - (107.0 - 114.0)).abs() < 1e-12);
        let dt = 1e-6;
        let step = fs_phs::step(&system, &initial, &[3e-4], dt).unwrap();
        let midpoint: Vec<_> = initial
            .iter()
            .zip(&step.x)
            .map(|(a, b)| f64::midpoint(*a, *b))
            .collect();
        let outward = flow.volume_flow_m3_s(&midpoint).unwrap();
        let volume_change: f64 = (0..cells)
            .map(|i| {
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                sign * (step.x[2 * i] - initial[2 * i])
            })
            .sum();
        let expected_volume = dt * (3e-4 - outward);
        assert!((volume_change - expected_volume).abs() < 1e-10 * expected_volume.abs());
        let travel_s: f64 = sections
            .iter()
            .zip(media)
            .map(|(s, m)| s.length / m.sound_speed)
            .sum();
        let omega = core::f64::consts::PI / (2.0 * travel_s);
        let resistance = 1.2 * omega.powi(2) / (2.0 * core::f64::consts::PI * 340.0);
        let expected_work = dt * resistance * outward.powi(2);
        assert!((step.dissipated - expected_work).abs() < 1e-10 * expected_work);
    }
}

#[test]
fn regional_foster_memory_uses_each_cells_own_transport_without_changing_its_neighbors() {
    let sections = [section(2, false); 2];
    let pin = ViscothermalPin {
        dynamic_viscosity: 1.8e-5,
        gamma: 1.4,
        prandtl: 0.71,
        foster_branches: 2,
    };
    let first = AcousticSectionMedium {
        density: 1.2,
        sound_speed: 340.0,
        viscothermal: Some(pin),
    };
    let mut media = [
        first,
        AcousticSectionMedium {
            density: 2.0,
            sound_speed: 280.0,
            viscothermal: Some(ViscothermalPin {
                foster_branches: 3,
                ..pin
            }),
        },
    ];
    let build = |media: &[AcousticSectionMedium]| {
        acoustic_chain_with_regional_outlet_flow(
            &sections,
            media,
            1.2,
            340.0,
            ChainOutlet::PressureRelease,
            1,
            None,
        )
        .unwrap()
    };
    let (before, _) = build(&media);
    media[1].viscothermal = Some(ViscothermalPin {
        dynamic_viscosity: 3.1e-5,
        gamma: 1.6,
        prandtl: 0.8,
        foster_branches: 3,
    });
    let (after, _) = build(&media);
    // Eight physical states plus 2*(2*2 + 2*3) local Foster states.
    assert_eq!(before.state_dim(), 28);
    assert_eq!(after.state_dim(), before.state_dim());
    let n = before.state_dim();
    let unaffected: Vec<_> = (0..4).chain(8..16).collect();
    let (_, before_r, _) = before.structure();
    let (_, after_r, _) = after.structure();
    for row in unaffected {
        assert_eq!(
            before_r[row * n..(row + 1) * n],
            after_r[row * n..(row + 1) * n]
        );
        let mut basis = vec![0.0; n];
        basis[row] = 1.0;
        assert_eq!(before.effort(&basis), after.effort(&basis));
    }
    assert_ne!(before_r, after_r);
    let step = fs_phs::step(&after, &vec![0.0; n], &[2e-4], 1e-5).unwrap();
    assert!(step.dissipated > 0.0);
    media[1].sound_speed = f64::NAN;
    assert!(
        acoustic_chain_with_regional_outlet_flow(
            &sections,
            &media,
            1.2,
            340.0,
            ChainOutlet::PressureRelease,
            1,
            None,
        )
        .is_err()
    );
    assert!(
        acoustic_chain_with_regional_outlet_flow(
            &sections,
            &[first],
            1.2,
            340.0,
            ChainOutlet::PressureRelease,
            1,
            None,
        )
        .is_err()
    );
}

#[test]
fn a_single_cone_uses_physical_volume_coordinates_while_the_legacy_source_keeps_its_psi_image() {
    for cells in [3, 4] {
        for outlet in [
            ChainOutlet::PressureRelease,
            ChainOutlet::Radiating(MouthFlange::Flanged),
        ] {
            let sections = [section(cells, true)];
            let original =
                acoustic_chain_outlet_wall(&sections, 1.2, 340.0, outlet, 1, &[], None, None)
                    .unwrap();
            assert_eq!(original.state_dim(), 2 * cells + 1);
            let step =
                fs_phs::step(&original, &vec![0.0; original.state_dim()], &[2e-4], 1e-5).unwrap();
            assert!(step.x.iter().any(|x| *x != 0.0));
            let (physical, flow) =
                acoustic_chain_with_outlet_flow(&sections, 1.2, 340.0, outlet, 1, &[], None, None)
                    .unwrap();
            assert_eq!(physical.state_dim(), 2 * cells);
            let step =
                fs_phs::step(&physical, &vec![0.0; physical.state_dim()], &[2e-4], 1e-5).unwrap();
            assert!(flow.volume_flow_m3_s(&step.x).unwrap() > 0.0);
        }
    }
}

#[test]
fn rigid_terminal_never_radiates_its_disconnected_flux_and_invalid_histories_refuse() {
    let (system, flow) = acoustic_chain_with_outlet_flow(
        &[section(3, true)],
        1.2,
        340.0,
        ChainOutlet::Rigid,
        1,
        &[],
        None,
        None,
    )
    .unwrap();
    let mut state = vec![1.0; system.state_dim()];
    assert_eq!(flow.volume_flow_m3_s(&state).unwrap(), 0.0);
    assert!(flow.volume_flow_m3_s(&state[1..]).is_err());
    state[0] = f64::NAN;
    assert!(flow.volume_flow_m3_s(&state).is_err());
}
