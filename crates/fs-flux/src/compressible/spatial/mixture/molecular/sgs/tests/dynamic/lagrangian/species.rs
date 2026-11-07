use super::*;
use fs_lbm::sgs3::dynamic::lagrangian::scalar::components::LagrangianDynamicComponents;

mod accepted;

fn config() -> MixtureLagrangianSpeciesSubgridConfig {
    MixtureLagrangianSpeciesSubgridConfig {
        momentum: memory_config().transport.model,
        heat: LagrangianScalarConfig {
            memory_time: 0.017,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        species: LagrangianComponentsConfig {
            memory_time: 0.023,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn mixing_field() -> MixtureEuler3d {
    let grid = CartesianGrid::new(9, 7, 5, 0.09, 0.14, 0.20).unwrap();
    field(material(0.0), grid, PERIODIC, |point| {
        let [x, y, _] = std::array::from_fn(|axis| {
            std::f64::consts::TAU * point[axis] / (grid.shape[axis] as f64 * grid.spacing[axis])
        });
        let fraction = 0.5 + 0.125 * (x + y).sin();
        (
            [fraction, 1.0 - fraction, 0.0],
            [0.3 * y.sin(), 0.0, 0.0],
            1.0 + 0.1 * (x + y).sin(),
        )
    })
}

fn independent_owner(formation: f64) -> (Coupled, ForcingWorkspace3d) {
    let flow = dynamic_field(formation);
    let molecular = operator(&flow, None);
    let work = flow.forcing_workspace(usize::MAX).unwrap();
    let bytes = Coupled::required_workspace_bytes(&flow, &molecular).unwrap()
        + Coupled::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow).unwrap()
        + work.workspace_bytes();
    let owner =
        Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, config(), bytes, || false)
            .unwrap();
    assert_eq!(owner.workspace_bytes() + work.workspace_bytes(), bytes);
    (owner, work)
}

fn scalar_fields(flow: &MixtureEuler3d) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::with_capacity(flow.cell_count());
    let mut velocity = Vec::with_capacity(flow.cell_count());
    let mut temperature = Vec::with_capacity(flow.cell_count());
    for cell in 0..flow.cell_count() {
        let primitive = flow
            .cell(cell)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap();
        rho.push(primitive.density);
        velocity.push(primitive.velocity);
        temperature.push(primitive.temperature);
    }
    (rho, velocity, temperature)
}

fn reference_components(flow: &MixtureEuler3d) -> LagrangianDynamicComponents {
    LagrangianDynamicComponents::new(
        flow.grid.shape,
        flow.grid.spacing,
        flow.gas().species_count(),
        usize::MAX,
    )
    .unwrap()
}

fn fits(
    owner: &Coupled,
) -> (
    Option<DynamicSgsReport>,
    Option<DynamicScalarReport>,
    Option<DynamicScalarReport>,
) {
    (
        owner.lagrangian_subgrid_report(),
        owner.lagrangian_subgrid_heat_report(),
        owner.lagrangian_subgrid_species_report(),
    )
}

#[test]
fn g1_lagrangian_heat_species_faces_survive_zero_momentum_fit() {
    let flow = mixing_field();
    let mut active = operator(&flow, None);
    let mut base = operator(&flow, None);
    active.sgs = Some(SubgridTransport::new_lagrangian_species(&flow, config()).unwrap());
    let bytes = active.workspace_bytes();
    let mut history = active
        .sgs
        .as_mut()
        .unwrap()
        .take_lagrangian_filter()
        .unwrap();
    assert_eq!(
        active.workspace_bytes(),
        bytes,
        "detached histories remain budgeted"
    );
    let width = history.filter_width();
    let trial = active
        .sgs
        .as_mut()
        .unwrap()
        .prepare_lagrangian(&mut history, &flow, &mut || false)
        .unwrap();
    let (momentum, heat, species) = trial.reports();
    assert_eq!(momentum.maximum_coefficient_squared, 0.0);
    assert!(heat.unwrap().maximum_coefficient > 0.0);
    assert!(species.unwrap().maximum_coefficient > 0.0);
    assert!(faces(&mut active, &flow).rate > faces(&mut base, &flow).rate);
    let ViscosityClosure::Lagrangian(model) = &active.sgs.as_ref().unwrap().closure else {
        unreachable!()
    };
    let mut heat_faces = 0;
    let mut species_faces = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let face = flow.grid.cell_face(axis, c);
            let heat_coefficient = model.face_heat_coefficient(left, right).unwrap().0;
            let species_coefficient = model.face_species_coefficient(left, right).unwrap().0;
            let rho_alpha = dynamic_face_viscosity_reference(
                &active,
                [heat_coefficient; 2],
                width,
                left,
                right,
                axis,
            );
            let rho_diffusion = dynamic_face_viscosity_reference(
                &active,
                [species_coefficient; 2],
                width,
                left,
                right,
                axis,
            );
            let pl = active.work.primitive[left];
            let pr = active.work.primitive[right];
            let distance = flow.grid.spacing[axis];
            let yl = flow.states()[7 * left] / pl.density;
            let yr = flow.states()[7 * right] / pr.density;
            let flux = -rho_diffusion * (yr - yl) / distance;
            let delta: [f64; 7] = std::array::from_fn(|k| {
                active.work.flux[axis][7 * face + k] - base.work.flux[axis][7 * face + k]
            });
            close(delta[0], flux, 4e-14);
            close(delta[1], -flux, 4e-14);
            assert_eq!(delta[2], 0.0);
            for momentum_flux in &delta[3..6] {
                assert_eq!(*momentum_flux, 0.0);
            }
            let temperature = Temperature::new(0.5 * pl.temperature + 0.5 * pr.temperature);
            let mut cp = 0.0;
            let mut enthalpy_flux = 0.0;
            for (k, source) in flow.gas().source().models().iter().take(2).enumerate() {
                let yl = flow.states()[7 * left + k] / pl.density;
                let yr = flow.states()[7 * right + k] / pr.density;
                cp += (0.5 * yl + 0.5 * yr)
                    * source
                        .evaluate_properties(temperature)
                        .unwrap()
                        .cp()
                        .value()
                    / source.molar_mass().value();
                let h = (source
                    .molar_internal_energy_change(
                        Temperature::new(flow.gas().domain().minimum_temperature),
                        temperature,
                    )
                    .unwrap()
                    .value()
                    + R * temperature.value())
                    / source.molar_mass().value();
                enthalpy_flux += h * delta[k];
            }
            let heat = -rho_alpha * cp * (pr.temperature - pl.temperature) / distance;
            close(delta[6], heat + enthalpy_flux, 8e-14);
            heat_faces += usize::from(heat.abs() > 1e-8);
            species_faces += usize::from(flux.abs() > 1e-8);
        }
    }
    assert!(heat_faces > 0 && species_faces > 0);
    drop(trial);
    active
        .sgs
        .as_mut()
        .unwrap()
        .restore_lagrangian_filter(history, None);
    assert_eq!(active.workspace_bytes(), bytes);
    let model = active.sgs.as_ref().unwrap();
    assert_eq!(model.lagrangian_sample_time(), None);
    assert_eq!(model.lagrangian_report(), None);
    assert_eq!(model.lagrangian_heat_report(), None);
    assert_eq!(model.lagrangian_species_report(), None);
}

#[test]
fn g0_lagrangian_species_storage_bounds_and_refused_history_are_atomic() {
    let flow = mixing_field();
    let molecular = operator(&flow, None);
    let additional =
        Coupled::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(additional, 448 * flow.cell_count());
    let total = Coupled::required_workspace_bytes(&flow, &molecular).unwrap() + additional;
    assert!(
        Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, config(), total - 1, || {
            false
        },)
        .is_err()
    );
    for limit in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            MixtureLagrangianSpeciesSubgridConfig {
                maximum_species_diffusivity_m2_s: limit,
                ..config()
            }
            .validate()
            .is_err()
        );
    }
    let mut forbidden = config();
    forbidden.species.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    forbidden = config();
    forbidden.heat.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let flow = mixing_field();
    let molecular = dynamic_operator(&flow, dynamic_config());
    assert!(
        Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, config(), usize::MAX, || {
            false
        },)
        .is_err()
    );

    let mut thermal_bound = config();
    thermal_bound.heat.maximum_diffusivity_to_molecular = 0.0;
    let mut species_bound = config();
    species_bound.maximum_species_diffusivity_m2_s = 0.0;
    let mut coefficient_bound = config();
    coefficient_bound.species.maximum_coefficient = 0.0;
    for cfg in [thermal_bound, species_bound, coefficient_bound] {
        let flow = mixing_field();
        let molecular = operator(&flow, None);
        let mut refused =
            Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, cfg, usize::MAX, || {
                false
            })
            .unwrap();
        let before = refused.flow().states().to_vec();
        let loads = refused.molecular().boundary_impulses.clone();
        let bytes = refused.workspace_bytes();
        for _ in 0..2 {
            assert!(
                refused
                    .step(Time::new(0.001), MolecularStepControl::default(), || false,)
                    .is_err()
            );
            assert_eq!(refused.flow().states(), before);
            assert_eq!(refused.flow().time(), Time::new(0.0));
            assert_eq!(refused.molecular().boundary_impulses, loads);
            assert_eq!(refused.last_step(), None);
            assert_eq!(refused.molecular().last_step(), None);
            assert_eq!(refused.lagrangian_subgrid_sample_time(), None);
            assert_eq!(fits(&refused), (None, None, None));
            assert_eq!(refused.workspace_bytes(), bytes);
        }
    }
}

#[test]
fn g0_lagrangian_species_uniform_fit_preserves_molecular_bits_and_clock_owner() {
    let uniform = || {
        field(
            material(0.0),
            CartesianGrid::new(4, 4, 4, 0.25, 0.25, 0.25).unwrap(),
            PERIODIC,
            |_| ([0.6, 0.4, 0.0], [0.1, -0.05, 0.025], 1.0),
        )
    };
    let flow = uniform();
    let molecular = operator(&flow, None);
    let mut base = Coupled::new(flow, molecular, usize::MAX, || false).unwrap();
    let flow = uniform();
    let molecular = operator(&flow, None);
    let mut active =
        Coupled::with_lagrangian_dynamic_species_sgs(flow, molecular, config(), usize::MAX, || {
            false
        })
        .unwrap();
    let bytes = active.workspace_bytes();
    for _ in 0..3 {
        let expected = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let actual = active
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            active
                .lagrangian_subgrid_report()
                .unwrap()
                .maximum_coefficient_squared,
            0.0
        );
        assert_eq!(
            active
                .lagrangian_subgrid_heat_report()
                .unwrap()
                .maximum_coefficient,
            0.0
        );
        assert_eq!(
            active
                .lagrangian_subgrid_species_report()
                .unwrap()
                .maximum_coefficient,
            0.0
        );
        for (actual, expected) in active.flow().states().iter().zip(base.flow().states()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
        assert_eq!(active.workspace_bytes(), bytes);
    }
    let (mut flow, mut molecular) = active.into_parts();
    let before = flow.states().to_vec();
    let time = flow.time();
    assert!(
        molecular
            .advance(
                &mut flow,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .is_err(),
        "the molecular pseudo-clock cannot independently age any history"
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time(), time);
}
