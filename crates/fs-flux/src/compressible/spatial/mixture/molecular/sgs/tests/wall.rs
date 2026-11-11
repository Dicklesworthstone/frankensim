use super::*;
use fs_lbm::conjugate::turbulence::law_of_the_wall_ratios;

fn wall_config() -> EquilibriumWallConfig {
    EquilibriumWallConfig {
        maximum_wall_distance_plus: 300.0,
        maximum_mach_number: 0.3,
    }
}

fn wall(velocity: [f64; 3]) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity,
        thermal: WallThermalCondition::Adiabatic,
    }
}

fn channel(normal: usize, velocity: [f64; 3]) -> MixtureEuler3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[normal] = WALL;
    field(
        material(0.0),
        CartesianGrid::new(4, 4, 4, 0.2, 0.3, 0.4).unwrap(),
        Boundaries3d {
            x: axes[0],
            y: axes[1],
            z: axes[2],
        },
        |_| ([0.4, 0.6, 0.0], velocity, 2.0),
    )
}

fn source_with_conductivity(flow: &MixtureEuler3d, conductivity: f64) -> SpatialTransport<ConstantTransport> {
    source_scaled_diffusion(flow, conductivity, 1.0)
}

/// `source_with_conductivity` with every binary diffusivity multiplied by `diffusion_scale`. The scale 1.0
/// returns the same bits (`x * 1.0` is exact), so every caller that passes it is unchanged; a smaller scale
/// lets a test see the viscous stability bound, which the species-diffusion bound (about 3e-2 here, against
/// the viscous 6.7e-4) otherwise exceeds on every face.
fn source_scaled_diffusion(
    flow: &MixtureEuler3d,
    conductivity: f64,
    diffusion_scale: f64,
) -> SpatialTransport<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "equilibrium-mixture-sgs-wall",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(1e-4),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        [0.0, 0.01, 0.012, 0.01, 0.0, 0.013, 0.012, 0.013, 0.0]
            .into_iter()
            .map(|diffusivity| diffusivity * diffusion_scale)
            .collect(),
    )
    .unwrap();
    SpatialTransport::new(flow.gas().clone(), law, || false).unwrap()
}

fn momentum() -> ChannelDynamicSgsConfig {
    ChannelDynamicSgsConfig {
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn heat() -> ChannelDynamicScalarConfig {
    ChannelDynamicScalarConfig {
        maximum_coefficient: 100.0,
        maximum_diffusivity_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn prepared(
    flow: &MixtureEuler3d,
    normal: usize,
    velocity: [f64; 3],
    kind: usize,
) -> MolecularTransport3d<ConstantTransport> {
    prepared_with_conductivity(flow, normal, velocity, kind, 0.0)
}

fn prepared_with_conductivity(
    flow: &MixtureEuler3d,
    normal: usize,
    velocity: [f64; 3],
    kind: usize,
    conductivity: f64,
) -> MolecularTransport3d<ConstantTransport> {
    prepared_scaled_diffusion(flow, normal, velocity, kind, conductivity, 1.0)
}

/// `prepared_with_conductivity` on a source whose binary diffusivities are scaled (`source_scaled_diffusion`).
fn prepared_scaled_diffusion(
    flow: &MixtureEuler3d,
    normal: usize,
    velocity: [f64; 3],
    kind: usize,
    conductivity: f64,
    diffusion_scale: f64,
) -> MolecularTransport3d<ConstantTransport> {
    let source = source_scaled_diffusion(flow, conductivity, diffusion_scale);
    match kind {
        0 => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_sgs(flow, &source).unwrap();
            MolecularTransport3d::with_sgs(
                flow,
                source,
                config(SgsModel::Smagorinsky { cs: 0.0 }),
                bytes,
                || false,
            )
            .unwrap()
            .with_wall(BoundarySide::ALL[2 * normal], wall(velocity))
            .unwrap()
            .with_wall(BoundarySide::ALL[2 * normal + 1], wall(velocity))
            .unwrap()
        }
        1 => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_dynamic_channel_thermal_sgs(
                    flow, &source,
                )
                .unwrap();
            MolecularTransport3d::with_dynamic_channel_thermal_sgs(
                flow,
                source,
                normal,
                [wall(velocity); 2],
                MixtureDynamicChannelThermalSubgridConfig {
                    momentum: momentum(),
                    heat: heat(),
                    turbulent_schmidt: 0.8,
                },
                bytes,
                || false,
            )
            .unwrap()
        }
        _ => {
            let bytes =
                MolecularTransport3d::required_workspace_bytes_with_channel_dynamic_species_sgs(
                    flow, &source,
                )
                .unwrap();
            MolecularTransport3d::with_channel_dynamic_species_sgs(
                flow,
                source,
                MixtureChannelDynamicSpeciesSubgridConfig {
                    transport: ChannelDynamicThermalSubgridConfig {
                        momentum: momentum(),
                        heat: heat(),
                    },
                    species: ChannelDynamicComponentsConfig {
                        maximum_coefficient: 100.0,
                        usage: SgsUse::Exploration,
                    },
                    maximum_species_diffusivity_m2_s: 10.0,
                },
                [wall(velocity); 2],
                bytes,
                || false,
            )
            .unwrap()
        }
    }
}

fn isothermal(velocity: [f64; 3], temperature: f64) -> MolecularWall3d {
    MolecularWall3d::NoSlip {
        velocity,
        thermal: WallThermalCondition::Isothermal(Temperature::new(temperature)),
    }
}

fn thermal_config(partial: [f64; 3], wall_temperature: f64) -> EquilibriumThermalWallConfig {
    // Independent manufactured NASA-9 formulas, at the ACTUAL wall/cell
    // midpoint temperature. The three species have different molecular masses.
    let t = 1.0 + 0.5 * wall_temperature;
    let cp = partial[0] * (3.5 + 0.2 * t) + 0.5 * partial[1] * (4.5 + 0.3 * t)
        + 2.0 * partial[2] * (4.5 + 0.2 * t);
    EquilibriumThermalWallConfig {
        momentum: wall_config(),
        heat: WallHeatConfig {
            turbulent_prandtl: cp * 1e-4, // k = 1 W/(m K)
            relative_tolerance: 1e-10,
            max_refinements: 14,
        },
    }
}

#[test]
fn equilibrium_sgs_thermal_wall_mixture_uses_face_calorics_composition_and_complete_bound() {
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[a] = 0.2;
        velocity[b] = -0.1;
        let wall_velocity = velocity.map(|u| 0.25 * u);
        for partial in [[0.4, 0.6, 0.0], [0.7, 0.3, 0.0]] {
            let template = channel(axis, velocity);
            let flow = field(material(0.0), template.grid.clone(), template.boundaries,
                |_| (partial, velocity, 2.0));
            let distance = 0.5 * flow.grid.spacing[axis];
            let ratio = law_of_the_wall_ratios(
                (0.15_f64.powi(2) + 0.075_f64.powi(2)).sqrt() * distance / 1e-4).1;
            let cv = partial[0] * (2.5 + 0.2 * 2.0) + 0.5 * partial[1] * (3.5 + 0.3 * 2.0);
            for kind in 0..3 {
                for wall_temperature in [1.8, 2.2] {
                    let config = thermal_config(partial, wall_temperature);
                    let mut op = prepared_with_conductivity(&flow, axis, wall_velocity, kind, 1.0);
                    let bytes = op.workspace_bytes();
                    for side in &BoundarySide::ALL[2 * axis..2 * axis + 2] {
                        op = op.with_wall(*side, isothermal(wall_velocity, wall_temperature)).unwrap()
                            .with_equilibrium_sgs_thermal_wall(*side, config).unwrap();
                        assert_eq!(op.equilibrium_sgs_thermal_wall(*side), Some(config));
                    }
                    let stage = faces(&mut op, &flow);
                    close(stage.rate, 4.0 * (ratio / cv)
                        * flow.grid.spacing.map(|h| h.recip().powi(2)).iter().sum::<f64>(), 3e-11);
                    assert_eq!(op.workspace_bytes(), bytes);
                    for upper in [false, true] {
                        let sign = if upper { 1.0 } else { -1.0 };
                        let mut face = [1; 3];
                        face[axis] = if upper { 4 } else { 0 };
                        let index = op.grid.cell_face(axis, face) * 7;
                        let flux = &op.work.flux[axis][index..index + 7];
                        assert_eq!(&flux[..3], &[0.0; 3]);
                        let work = flux[3 + a] * wall_velocity[a] + flux[3 + b] * wall_velocity[b];
                        close(flux[6], work - ratio * sign * (wall_temperature - 2.0) / distance, 3e-11);
                    }
                }
            }
        }
    }
}

#[test]
fn equilibrium_sgs_thermal_wall_mixture_interval_conserves_inventory_and_boundary_heat() {
    for wall_temperature in [1.8, 2.2] {
        let mut flow = channel(1, [0.4, 0.0, 0.0]);
        let mut op = prepared_with_conductivity(&flow, 1, [0.1, 0.0, 0.0], 0, 1.0);
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            op = op.with_wall(side, isothermal([0.1, 0.0, 0.0], wall_temperature)).unwrap()
                .with_equilibrium_sgs_thermal_wall(side, thermal_config([0.4, 0.6, 0.0], wall_temperature)).unwrap();
        }
        let before = flow.totals();
        let report = op.advance(&mut flow, Time::new(0.001), DissipationControl::default(), || false).unwrap();
        let heat = report.boundary_change.energy - 0.1 * report.boundary_change.momentum[0];
        assert_eq!(heat.signum(), (wall_temperature - 2.0).signum());
        close(flow.totals().energy - before.energy, report.boundary_change.energy, 3e-12);
        for species in 0..3 { close(flow.totals().species_masses[species], before.species_masses[species], 3e-12); }
        assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
        balance(&report);
        let states = flow.states().to_vec();
        let loads = op.boundary_impulses.clone();
        let report = op.last_step().cloned();
        // The lower wall has already prepared valid flux when the upper wall
        // exhausts its explicit quadrature work budget. Publish none of it.
        op = op.with_equilibrium_sgs_thermal_wall(BoundarySide::YMax,
            EquilibriumThermalWallConfig { heat: WallHeatConfig {
                turbulent_prandtl: 0.9, relative_tolerance: 1e-12, max_refinements: 3,
            }, ..thermal_config([0.4, 0.6, 0.0], wall_temperature) }).unwrap();
        assert!(op.advance(&mut flow, Time::new(0.001), DissipationControl::default(), || false).is_err());
        assert_eq!(flow.states(), states);
        assert_eq!(op.boundary_impulses, loads);
        assert_eq!(op.last_step(), report.as_ref());
    }
}

#[test]
fn authored_equilibrium_wall_keeps_its_complete_native_traction_operator() {
    let velocity = [0.0, 0.02, 0.0];
    let flow = channel(0, velocity);
    let original = prepared(&flow, 0, [0.0; 3], 0);
    assert!(original.uses_compatible_molecular_stress());
    let modeled = original
        .with_equilibrium_sgs_wall(BoundarySide::XMin, wall_config())
        .unwrap();
    assert!(!modeled.uses_compatible_molecular_stress());
}

#[test]
fn equilibrium_sgs_wall_mixture_faces_preserve_species_and_native_work_for_local_and_dynamic_models()
 {
    for axis in 0..3 {
        let [a, b] = tangents(axis);
        let mut velocity = [0.0; 3];
        velocity[a] = 0.2;
        velocity[b] = -0.1;
        let wall_velocity = velocity.map(|u| 0.25 * u);
        let flow = channel(axis, velocity);
        let distance = 0.5 * flow.grid.spacing[axis];
        // The fixture's dynamic viscosity (`source_scaled_diffusion`; its bulk viscosity is 0) and density (`channel`).
        let mu = 1e-4;
        let rho = 0.4 + 0.6;
        let (eddy_ratio, ratio) =
            law_of_the_wall_ratios((0.15_f64.powi(2) + 0.075_f64.powi(2)).sqrt() * distance / mu);
        // The claim is that attaching the wall adds its nonlinear Jacobian term, 4 (tangent - molecular) / rho, to
        // the wall face's diffusivity and so to the stability rate. The rate is the
        // maximum over faces; with the fixture's own binary diffusivities (0.01 to 0.013) the species-diffusion
        // bound (about 3.4e-2) is that maximum on every face, 30 times the wall face's 1.1e-3, and no wall term
        // can show in it. With the diffusivities scaled by 1e-3 the viscous bound is the maximum.
        let diffusion_scale = 1e-3;
        let scaling = 4.0 * flow.grid.spacing.iter().map(|h| h.recip().powi(2)).sum::<f64>();
        let mut base =
            prepared_scaled_diffusion(&flow, axis, wall_velocity, 0, 0.0, diffusion_scale);
        let base_rate = faces(&mut base, &flow).rate;
        // The exact identity below holds only if the wall face's diffusivity before the Jacobian term is the base
        // maximum, the viscous bound (4/3 mu + 2 |bulk - 2/3 mu| + 4 mu) / rho with bulk 0, which every face has.
        // Roundings: three products and two sums for the bound, a division, and the scaling (three squares, two
        // sums, a product) times it: 16 eps relative bounds them.
        let floor = (4.0 / 3.0 * mu + 2.0 * (2.0 / 3.0 * mu) + 4.0 * mu) / rho * scaling;
        assert!(
            (base_rate - floor).abs() <= 16.0 * f64::EPSILON * floor,
            "base rate {base_rate:e} is not the viscous bound {floor:e}: species diffusion still binds"
        );
        for kind in 0..3 {
            let mut op = prepared_scaled_diffusion(
                &flow,
                axis,
                wall_velocity,
                kind,
                0.0,
                diffusion_scale,
            );
            let bytes = op.workspace_bytes();
            for upper in [false, true] {
                let side = BoundarySide::ALL[2 * axis + usize::from(upper)];
                op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
                assert_eq!(op.equilibrium_sgs_wall(side), Some(wall_config()));
            }
            let walled_rate = faces(&mut op, &flow).rate;
            assert!(
                walled_rate > base_rate,
                "walled rate {walled_rate:e} does not exceed the base rate {base_rate:e}"
            );
            if kind == 0 {
                // Smagorinsky with cs = 0 adds no eddy viscosity, so the wall face is the maximum and the walled
                // rate is the base rate plus the wall term times the scaling. The wall term's tangent is the
                // Jacobian of the discrete wall flux, d(secant Re)/dRe, not the interior eddy ratio: with
                // sigma = secant, s = dy+/du+ = 1 + eddy_ratio, Re = sigma u+^2 and dRe/du+ = Re/u+ + u+ s, it is
                // 2 sigma s / (sigma + s) (the helper's two public outputs; the operator's own tangent comes from a
                // differently ordered Newton and agrees with it to about 3e-16 relative). Roundings: that is a
                // sum, two products and a division; the `- 1` cancels about half of it (a factor 1.85 on its
                // relative error), and the term is 0.41 of the rate; the operator's (bound + term) * scaling is
                // one sum and one product, this side's base + term * scaling two products and a sum: 16 eps
                // relative bounds them. If the wall term were dropped, the walled rate would equal the base rate
                // and both assertions would fail.
                let jacobian = 2.0 * ratio * (1.0 + eddy_ratio) / (ratio + 1.0 + eddy_ratio);
                let expected = base_rate + 4.0 * mu * (jacobian - 1.0) / rho * scaling;
                assert!(
                    (walled_rate - expected).abs() <= 16.0 * f64::EPSILON * expected,
                    "walled rate {walled_rate:e} against base + wall term {expected:e}"
                );
            }
            assert_eq!(op.workspace_bytes(), bytes);
            for upper in [false, true] {
                let sign = if upper { 1.0 } else { -1.0 };
                let mut face = [1; 3];
                face[axis] = if upper { 4 } else { 0 };
                let index = op.grid.cell_face(axis, face) * 7;
                let flux = &op.work.flux[axis][index..index + 7];
                assert_eq!(&flux[..3], &[0.0; 3]);
                close(flux[3 + axis], 0.0, 1e-14);
                for t in [a, b] {
                    let momentum =
                        1e-4 * ratio * sign * (velocity[t] - wall_velocity[t]) / distance;
                    close(flux[3 + t], momentum, 2e-14);
                }
                close(
                    flux[6],
                    flux[3 + a] * wall_velocity[a] + flux[3 + b] * wall_velocity[b],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn equilibrium_sgs_wall_mixture_braking_conserves_species_and_total_energy_with_wall_work() {
    for speed in [0.0, 0.1] {
        let mut base = channel(1, [0.4, 0.0, 0.0]);
        let mut active = channel(1, [0.4, 0.0, 0.0]);
        let mut base_op = prepared(&base, 1, [speed, 0.0, 0.0], 0);
        let mut op = prepared(&active, 1, [speed, 0.0, 0.0], 0);
        for side in [BoundarySide::YMin, BoundarySide::YMax] {
            op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
        }
        let before = active.totals();
        base_op
            .advance(
                &mut base,
                Time::new(0.02),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let report = op
            .advance(
                &mut active,
                Time::new(0.02),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        assert!(active.totals().momentum[0] < base.totals().momentum[0]);
        close(
            report.boundary_change.energy,
            speed * report.boundary_change.momentum[0],
            1e-13,
        );
        close(
            active.totals().energy - before.energy,
            report.boundary_change.energy,
            2e-13,
        );
        for i in 0..3 {
            close(
                active.totals().species_masses[i],
                before.species_masses[i],
                1e-13,
            );
        }
        assert!(active.states().chunks_exact(7).all(|row| row[2] == 0.0));
        balance(&report);
    }
}

#[test]
fn equilibrium_sgs_wall_mixture_refusal_preserves_species_loads_and_report() {
    let mut flow = channel(1, [0.4, 0.0, 0.0]);
    let mut op = prepared(&flow, 1, [0.0; 3], 0);
    for side in [BoundarySide::YMin, BoundarySide::YMax] {
        op = op.with_equilibrium_sgs_wall(side, wall_config()).unwrap();
    }
    op.advance(
        &mut flow,
        Time::new(0.01),
        DissipationControl::default(),
        || false,
    )
    .unwrap();
    let states = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    let report = op.last_step().cloned();
    op = op
        .with_equilibrium_sgs_wall(
            BoundarySide::YMax,
            EquilibriumWallConfig {
                maximum_mach_number: 1e-6,
                ..wall_config()
            },
        )
        .unwrap();
    assert!(
        op.advance(
            &mut flow,
            Time::new(0.01),
            DissipationControl::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), states);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), report.as_ref());
    assert!(
        op.with_wall(
            BoundarySide::YMin,
            MolecularWall3d::NoSlip {
                velocity: [0.0; 3],
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
            }
        )
        .is_err()
    );
}
