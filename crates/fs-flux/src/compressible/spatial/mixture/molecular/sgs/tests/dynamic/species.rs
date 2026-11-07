use super::*;

fn config() -> MixtureDynamicSpeciesSubgridConfig {
    MixtureDynamicSpeciesSubgridConfig {
        transport: DynamicThermalSubgridConfig {
            momentum: dynamic_config().transport.model,
            heat: DynamicScalarConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient: 100.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        },
        species: DynamicComponentsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn species_field(formation: f64) -> MixtureEuler3d {
    let grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.20).unwrap();
    field(material(formation), grid, PERIODIC, |point| {
        let [x, y, _] = std::array::from_fn(|axis| {
            std::f64::consts::TAU * point[axis] / (grid.shape[axis] as f64 * grid.spacing[axis])
        });
        let fraction = 0.5 + 0.125 * (x + y).sin();
        (
            [fraction, 1.0 - fraction, 0.0],
            [0.3 * y.sin(), 0.0, 0.0],
            1.0 + 0.1 * (x + y).cos(),
        )
    })
}

fn prepared(
    flow: &MixtureEuler3d,
    cfg: MixtureDynamicSpeciesSubgridConfig,
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_species_sgs(flow, &source)
            .unwrap();
    let base = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
    assert_eq!(bytes - base, 352 * flow.cell_count());
    let op =
        MolecularTransport3d::with_dynamic_species_sgs(flow, source, cfg, bytes, || false).unwrap();
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(op.dynamic_species_subgrid_config(), Some(cfg));
    assert_eq!(op.dynamic_thermal_subgrid_config(), None);
    assert_eq!(op.dynamic_subgrid_config(), None);
    assert_eq!(op.subgrid_config(), None);
    op
}

fn fixed_schmidt(flow: &MixtureEuler3d) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_thermal_sgs(flow, &source)
            .unwrap();
    MolecularTransport3d::with_dynamic_thermal_sgs(
        flow,
        source,
        MixtureDynamicThermalSubgridConfig {
            transport: config().transport,
            turbulent_schmidt: 0.8,
        },
        bytes,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_dynamic_species_flux_and_source_enthalpy_survive_zero_eddy_viscosity() {
    let flow = species_field(0.0);
    let mut active = prepared(&flow, config());
    let mut fixed = fixed_schmidt(&flow);
    let rate = faces(&mut active, &flow).rate;
    let base_rate = faces(&mut fixed, &flow).rate;
    assert!(rate > base_rate);
    assert!(
        active
            .sgs
            .as_ref()
            .unwrap()
            .dynamic_coefficient_squared
            .iter()
            .all(|&mu| mu == 0.0)
    );
    let ViscosityClosure::DynamicSpecies(species) = &active.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    let mut nonzero = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let face = flow.grid.cell_face(axis, c);
            let pl = active.work.primitive[left];
            let pr = active.work.primitive[right];
            let mass_diffusion = dynamic_face_viscosity_reference(
                &active,
                [species.face_species_coefficient(left, right); 2],
                species.filter_width(),
                left,
                right,
                axis,
            );
            let gradient = (flow.states()[7 * right] / pr.density
                - flow.states()[7 * left] / pl.density)
                / flow.grid.spacing[axis];
            let flux = -mass_diffusion * gradient;
            let delta: [f64; 7] = std::array::from_fn(|k| {
                active.work.flux[axis][7 * face + k] - fixed.work.flux[axis][7 * face + k]
            });
            close(delta[0], flux, 4e-14);
            close(delta[1], -flux, 4e-14);
            assert_eq!(delta[2], 0.0);
            for k in 3..6 {
                assert_eq!(delta[k], 0.0);
            }
            let at = Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature);
            let h: Vec<_> = flow
                .gas()
                .source()
                .models()
                .iter()
                .take(2)
                .map(|model| {
                    (model
                        .molar_internal_energy_change(
                            Temperature::new(flow.gas().domain().minimum_temperature),
                            at,
                        )
                        .unwrap()
                        .value()
                        + R * at.value())
                        / model.molar_mass().value()
                })
                .collect();
            close(delta[6], (h[0] - h[1]) * flux, 4e-14);
            nonzero += usize::from(flux.abs() > 1e-8);
        }
    }
    assert!(nonzero > 0);
}

#[test]
fn g1_independent_heat_and_species_survive_centered_null_face_shear() {
    let mut flow = centered_null_field(true);
    let mut cfg = config();
    cfg.transport.momentum.averaging = DynamicAveraging::Domain;
    cfg.transport.heat.averaging = DynamicAveraging::Domain;
    cfg.species.averaging = DynamicAveraging::Domain;
    let mut active = prepared(&flow, cfg);
    let mut base = operator(&flow, None);
    let active_rate = faces(&mut active, &flow).rate;
    let base_rate = faces(&mut base, &flow).rate;
    assert!(active_rate > base_rate);
    let left = flow.grid.index([1, 1, 1]);
    let right = flow.grid.index([2, 1, 1]);
    let sgs = active.sgs.as_ref().unwrap();
    let ViscosityClosure::DynamicSpecies(model) = &sgs.closure else {
        unreachable!()
    };
    for cell in [left, right] {
        assert!(active.work.gradient[cell].iter().flatten().all(|g| g.abs() < 1e-15));
        assert!(model.face_heat_coefficient(cell, cell).0 > 0.0);
        assert!(model.face_species_coefficient(cell, cell) > 0.0);
    }
    let pl = active.work.primitive[left];
    let pr = active.work.primitive[right];
    let coefficient = model.face_species_coefficient(left, right);
    // At this face du_y/dx=-0.02, so |S|=0.02 despite both centered
    // tensors vanishing. Density and coefficients are averaged separately.
    let rho = 0.5 * pl.density + 0.5 * pr.density;
    let mass_diffusion = rho * coefficient * 0.02;
    let j = -mass_diffusion
        * (flow.states()[7 * right] / pr.density - flow.states()[7 * left] / pl.density);
    let face = flow.grid.cell_face(0, [2, 1, 1]) * 7;
    let delta: [f64; 7] = std::array::from_fn(|k| {
        active.work.flux[0][face + k] - base.work.flux[0][face + k]
    });
    assert!(j.abs() > 1e-8);
    close(delta[0], j, 3e-13);
    close(delta[1], -j, 3e-13);
    assert_eq!(delta[2], 0.0);
    let t = 0.5 * pl.temperature + 0.5 * pr.temperature;
    let mut cp = 0.0;
    let mut species_energy = 0.0;
    for (k, source) in flow.gas().source().models().iter().take(2).enumerate() {
        let yl = flow.states()[7 * left + k] / pl.density;
        let yr = flow.states()[7 * right + k] / pr.density;
        cp += 0.5 * (yl + yr) / source.molar_mass().value()
            * source.evaluate_properties(Temperature::new(t)).unwrap().cp().value();
        let h = (source.molar_internal_energy_change(
            Temperature::new(flow.gas().domain().minimum_temperature), Temperature::new(t),
        ).unwrap().value() + R * t) / source.molar_mass().value();
        species_energy += h * delta[k];
    }
    let heat = -rho * model.face_heat_coefficient(left, right).0 * 0.02 * cp
        * (pr.temperature - pl.temperature);
    assert!(heat.abs() > 1e-8);
    let stress_work: f64 = (0..3).map(|k| {
        (0.5 * pl.velocity[k] + 0.5 * pr.velocity[k]) * delta[3 + k]
    }).sum();
    close(delta[6], heat + species_energy + stress_work, 3e-13);
    let bytes = active.workspace_bytes() + flow.workspace_bytes();
    let report = active.advance(
        &mut flow, Time::new(0.001), DissipationControl::default(), || false,
    ).unwrap();
    balance(&report);
    assert_eq!(active.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
}

#[test]
fn g4_species_face_bound_refuses_when_centered_diffusivity_is_admitted() {
    use fs_lbm::sgs3::dynamic::scalar::PeriodicDynamicScalar;
    use fs_lbm::sgs3::dynamic::scalar::components::PeriodicDynamicComponents;

    let mut flow = centered_null_field(true);
    let mut cfg = config();
    cfg.transport.momentum.averaging = DynamicAveraging::Domain;
    cfg.transport.heat.averaging = DynamicAveraging::Domain;
    cfg.species.averaging = DynamicAveraging::Domain;
    let mut active = prepared(&flow, cfg);
    let controls = DissipationControl::default();
    let accepted = active.advance(&mut flow, Time::new(0.001), controls, || false).unwrap();
    faces(&mut active, &flow);
    let rho: Vec<_> = active.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = active.work.primitive.iter().map(|p| p.velocity).collect();
    let temperature: Vec<_> = active.work.primitive.iter().map(|p| p.temperature).collect();
    let mut core = PeriodicDynamicSgs::new_cartesian(
        flow.grid.shape, flow.grid.spacing, usize::MAX,
    ).unwrap();
    let mut heat = PeriodicDynamicScalar::new(flow.grid.shape, usize::MAX).unwrap();
    let mut species = PeriodicDynamicComponents::new(flow.grid.shape, 3, usize::MAX).unwrap();
    let (_, _, centered) = core.prepare_with_scalar_and_components(
        cfg.transport.momentum,
        cfg.transport.heat,
        cfg.species,
        &rho,
        &velocity,
        &temperature,
        &mut heat,
        |k, cell| flow.states()[7 * cell + k] / rho[cell],
        &mut species,
        || false,
    ).unwrap();
    let mut maximum_face: f64 = 0.0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let weight = dynamic_face_viscosity_reference(
                &active,
                [species.coefficient(left).unwrap(), species.coefficient(right).unwrap()],
                core.filter_width(),
                left,
                right,
                axis,
            );
            maximum_face = maximum_face.max(weight / (0.5 * rho[left] + 0.5 * rho[right]));
        }
    }
    assert!(maximum_face > 1.1 * centered.maximum_diffusivity);
    let bound = 0.5 * centered.maximum_diffusivity + 0.5 * maximum_face;
    assert!(centered.maximum_diffusivity < bound && bound < maximum_face);
    let ViscosityClosure::DynamicSpecies(model) = &mut active.sgs.as_mut().unwrap().closure else {
        unreachable!()
    };
    model.config.maximum_species_diffusivity_m2_s = bound;
    let state = flow.states().to_vec();
    let loads = active.boundary_impulses.clone();
    assert_eq!(
        active.advance(&mut flow, Time::new(0.001), controls, || false),
        Err(DissipationError::Invalid(
            "dynamic species face diffusivity exceeds its physical admission bound",
        )),
    );
    assert_eq!(flow.states(), state);
    assert_eq!(active.boundary_impulses, loads);
    assert_eq!(active.last_step(), Some(&accepted));
}

#[test]
fn g3_dynamic_species_updates_conserve_mix_and_preserve_source_gauge() {
    let mut flow = species_field(0.0);
    let mut shifted = species_field(1e12);
    let mut base = species_field(0.0);
    let mut active = prepared(&flow, config());
    let mut shifted_op = prepared(&shifted, config());
    let mut fixed = fixed_schmidt(&base);
    let bytes = active.workspace_bytes() + flow.workspace_bytes();
    let controls = DissipationControl::default();
    let report = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    balance(&report);
    shifted_op
        .advance(&mut shifted, Time::new(0.001), controls, || false)
        .unwrap();
    fixed
        .advance(&mut base, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(flow.states(), shifted.states());
    let variance = |flow: &MixtureEuler3d| {
        flow.states()
            .chunks_exact(7)
            .map(|row| {
                let rho = sum(row[..3].iter().copied());
                rho * (row[0] / rho - 0.5).powi(2)
            })
            .sum::<f64>()
    };
    assert!(variance(&flow) < variance(&base));
    assert!(flow.states().chunks_exact(7).all(|row| row[2] == 0.0));
    assert_eq!(active.workspace_bytes() + flow.workspace_bytes(), bytes);
}

#[test]
fn g0_dynamic_species_zero_fit_clears_cache_and_preserves_molecular_bits() {
    let flow = species_field(0.0);
    let mut active = prepared(&flow, config());
    faces(&mut active, &flow);
    let uniform = field(material(0.0), flow.grid, PERIODIC, |_| {
        ([0.5, 0.5, 0.0], [0.0; 3], 1.0)
    });
    let mut base = operator(&uniform, None);
    faces(&mut active, &uniform);
    faces(&mut base, &uniform);
    let ViscosityClosure::DynamicSpecies(species) = &active.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    assert!((0..flow.cell_count()).all(|cell| species.face_species_coefficient(cell, cell) == 0.0));
    assert_eq!(active.work.flux, base.work.flux);
}

#[test]
fn g4_dynamic_species_bounds_and_late_cancellation_preserve_accepted_state() {
    let mut flow = species_field(0.0);
    let source = original(&flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_species_sgs(&flow, &source)
            .unwrap();
    assert!(
        MolecularTransport3d::with_dynamic_species_sgs(&flow, source, config(), bytes - 1, || {
            false
        },)
        .is_err()
    );
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut cfg = config();
        cfg.maximum_species_diffusivity_m2_s = value;
        assert!(cfg.validate().is_err());
    }
    let mut forbidden = config();
    forbidden.species.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let mut cfg = config();
    cfg.maximum_species_diffusivity_m2_s = 0.0;
    let mut active = prepared(&flow, cfg);
    let before = flow.states().to_vec();
    let controls = DissipationControl::default();
    assert!(
        active
            .advance(&mut flow, Time::new(0.001), controls, || false)
            .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.last_step(), None);
    let ViscosityClosure::DynamicSpecies(species) = &mut active.sgs.as_mut().unwrap().closure
    else {
        unreachable!()
    };
    species.config = config();

    let mut reference_flow = species_field(0.0);
    let mut reference = prepared(&reference_flow, config());
    let mut polls = 0;
    let report = reference
        .advance(&mut reference_flow, Time::new(0.001), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    let loads = active.boundary_impulses.clone();
    let mut cancelled_polls = 0;
    assert_eq!(
        active.advance(&mut flow, Time::new(0.001), controls, || {
            cancelled_polls += 1;
            cancelled_polls == polls - 1
        }),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(active.boundary_impulses, loads);
    assert_eq!(active.last_step(), None);
    let retry = active
        .advance(&mut flow, Time::new(0.001), controls, || false)
        .unwrap();
    assert_eq!(retry, report);
    assert_eq!(flow.states(), reference_flow.states());
}
