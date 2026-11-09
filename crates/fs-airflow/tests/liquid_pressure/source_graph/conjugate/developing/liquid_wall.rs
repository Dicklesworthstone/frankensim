//! Same liquid source at the actual solved conducting-wall temperature.
use super::*;

fn liquid_h(rate: f64, surface: f64) -> f64 {
    let (_, diameter, _) = dimensions(false);
    let re: f64 = 20_000.0;
    let pr: f64 = 8.0;
    let friction = (1.8 * re.log10() - 1.5).powi(-2);
    let developed = friction / 8.0 * (re - 1000.0) * pr
        / (1.0 + 12.7 * (friction / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
    let start = if rate > 0.0 { 0.0 } else { 0.5 };
    let a = (start / diameter).cbrt();
    let b = ((start + 0.5) / diameter).cbrt();
    let wall_pr = 2000.0 * (0.005 - 1e-5 * surface) / (0.005 * surface - 1.0);
    developed * (1.0 + 1.0 / (a * a + a * b + b * b))
        * (pr / wall_pr).powf(0.11) * 0.5 / diameter
}

fn independent_surface(rate: f64, wall: f64) -> f64 {
    let (mut lo, mut hi) = (wall.min(300.0), wall.max(300.0));
    for _ in 0..90 {
        let surface = 0.5 * (lo + hi);
        let residual = 0.001 / DISTANCE * (wall - surface) * (wall + surface)
            - liquid_h(rate, surface) * (surface - 300.0);
        if residual > 0.0 { lo = surface; } else { hi = surface; }
    }
    0.5 * (lo + hi)
}

#[test]
fn g1_g3_liquid_wall_correction_uses_solved_surface_and_one_reciprocal_native_heat() {
    context(|cx| {
        let card = fluid_card();
        let (law, caloric) = source_pair(&card);
        let solid = solid(true);
        let wall = wall_caloric(&solid);
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        let layer = layer(&solid, &wall, &point, [200.0, 500.0]);
        let (_, diameter, area) = dimensions(false);
        let magnitude = 20_000.0 * 0.002 * PIPE_AREA / (1000.0 * diameter);
        for stored_wall in [400.0, 280.0, 300.0] {
            let initial = fixture(&caloric, &wall, stored_wall, area);
            let graph = graph();
            let inventory = bind_inventory(cx, &graph, &initial);
            let laws = laws(1e6);
            let owner = LiquidGraphSources::try_new(
                cx, &law, &caloric, &initial, SourceGraphBinding {
                    graph: &graph, inventory: &inventory, laws: &laws,
                    acceleration_m_s2: [0.0; 3],
                }, caps(),
            ).unwrap();
            let mut declared = binding(&layer, false, ConvectionPropertyTemperature::Reservoir);
            declared.section.duct = DevelopingDuct::TurbulentCircularLiquidWall {
                diameter: Length::new(diameter),
            };
            let prepared = PreparedGraphDevelopingConjugateLinks::new(
                &initial, vec![declared], 1, || Ok(()),
            ).unwrap();
            for sign in [-1.0, 1.0] {
                let q = sign * magnitude;
                let flows = [q, 0.0, q];
                let source = owner.resolve(cx, &initial, &flows).unwrap();
                let pressure = A * 0.002 * q;
                let report = native_report(
                    cx, &source.graph, &inventory, flows,
                    [pressure, 0.0, pressure], flows, 0.01,
                );
                let heat = owner.evaluate_developing_conjugate_links(
                    cx, &initial, &source.sources, &report, &prepared, config(),
                ).unwrap();
                let row = &heat.samples()[0];
                let film = &row.layer.film;
                let surface = independent_surface(q, stored_wall);
                let h = liquid_h(q, surface);
                let power = area * h * (surface - 300.0);
                near(row.layer.surface_temperature_k, surface);
                near(film.wall_temperature_k, surface);
                near(film.surface.as_ref().unwrap().thermal.temperature_k(), surface);
                near(row.wall_temperature_k, stored_wall);
                near(film.fluid_temperature_k, 300.0);
                near(film.property_temperature_k, 300.0);
                near(film.section.downstream().groups()["Re"], 20_000.0);
                near(film.section.downstream().groups()["Pr"], 8.0);
                near(film.coefficient.value.value(), h);
                near(film.conductance_w_k, area * h);
                assert_eq!(film.section.upstream().is_some(), sign < 0.0);
                let wall_source = film.surface.as_ref().unwrap();
                assert_eq!(wall_source.material.material().card_identity(), card.content_hash());
                near(wall_source.material.material().query().axes()["pressure"], 100000.0);
                let pr_wall = wall_source.thermal.transport().prandtl();
                near(film.coefficient.sensitivity.d_qoi["Pr_wall"], -0.11 * h / pr_wall);
                assert!((row.applied_heat_into_fluid_w - power).abs() <= 3e-10);
                let conducted = area / DISTANCE * 0.001
                    * (stored_wall - surface) * (stored_wall + surface);
                assert!((row.applied_heat_into_fluid_w - conducted).abs() <= 3e-10);
                assert!(row.layer.heat_residual_w.abs() <= 3e-10);
                assert!(Arc::ptr_eq(&row.source_line, layer.source_line()));
                if stored_wall == 300.0 {
                    assert_eq!(row.applied_heat_into_fluid_w, 0.0);
                    near(row.layer.conductance_w_k,
                        1.0 / (1.0 / (area * h) + DISTANCE / (area * 0.6)));
                } else {
                    assert!((surface - stored_wall).abs() > 1.0);
                    assert!((h / liquid_h(q, stored_wall) - 1.0).abs() > 0.01,
                        "using the stored wall as the property surface must change this fixture");
                }
                let mut links = BTreeMap::new();
                heat.apply_to(&mut links, || Ok(())).unwrap();
                let trial = BackwardEulerStep::new(&initial, 0.01, || Ok(())).unwrap()
                    .evaluate(&[0.0; 7], |_| Ok::<_, DepositionError>(StepCoefficients {
                        sources_w: BTreeMap::new(), conductances_w_k: links, payload: (),
                    }), || Ok(())).unwrap();
                near(trial.residual_j()[2], -0.01 * row.applied_heat_into_fluid_w);
                near(trial.residual_j()[6], 0.01 * row.applied_heat_into_fluid_w);
                near(trial.residual_j().iter().sum(), 0.0);
            }
            assert_eq!(initial.time_s(), 0.0);
            assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
        }
    });
}
