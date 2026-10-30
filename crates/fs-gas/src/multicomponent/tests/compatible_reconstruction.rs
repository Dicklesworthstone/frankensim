use super::*;

#[test]
fn g3_density_shocks_preserve_uniform_inert_and_linear_mass_fractions_on_both_meshes() {
    for nonuniform in [false, true] {
        let gas = gas(["A", "B", "C", "D", "inert"]
            .into_iter()
            .map(|name| source(name, R, 3.5, 0.0, 0.0))
            .collect());
        let count = 80;
        let positions: Vec<_> = (0..=count)
            .map(|i| {
                let x = i as f64 / count as f64;
                if nonuniform { x.powf(1.3) } else { x }
            })
            .collect();
        let area = 0.7;
        let volumes: Vec<_> = positions.windows(2).map(|x| area * (x[1] - x[0])).collect();
        let mesh = Mesh1d::new(positions.clone(), vec![area; count + 1], volumes).unwrap();
        let states = positions
            .windows(2)
            .map(|pair| {
                let x = 0.5 * (pair[0] + pair[1]);
                let phase = 2.0 * std::f64::consts::PI * x;
                let a = 0.20 + 0.08 * phase.sin();
                let c = 0.10 + 0.03 * (2.0 * phase + 0.4).sin();
                // This independent linear combination stays constant under
                // advection: YA + 2YB + 3YC = 1.1, while Y_inert=.3.
                let b = 0.55 - 0.5 * a - 1.5 * c;
                let d = 0.7 - a - b - c;
                let density = if x < 0.5 { 2.0 } else { 0.4 };
                let pressure = if x < 0.5 { 2.5 } else { 0.4 };
                let partial = [a, b, c, d, 0.3].map(|y| density * y);
                gas.conserved(&partial, 0.6, pressure).unwrap()
            })
            .collect();
        let mut flow =
            MixtureEuler1d::new_on_mesh(gas, mesh, states, Boundaries::Periodic, || false).unwrap();
        let initial = flow.totals();
        run(&mut flow, 0.08);
        let mut inert_error: f64 = 0.0;
        let mut linear_error: f64 = 0.0;
        for i in 0..flow.cell_count() {
            let q = flow.cell(i).unwrap();
            let rho = q.density();
            inert_error = inert_error.max((q.partial_densities[4] / rho - 0.3).abs());
            let linear = (q.partial_densities[0]
                + 2.0 * q.partial_densities[1]
                + 3.0 * q.partial_densities[2])
                / rho;
            linear_error = linear_error.max((linear - 1.1).abs());
        }
        eprintln!(
            "compatible species reconstruction nonuniform={nonuniform}: max inert error={inert_error:e}, linear error={linear_error:e}"
        );
        assert!(inert_error < 2e-13, "inert fraction drift {inert_error:e}");
        assert!(
            linear_error < 5e-13,
            "linear composition drift {linear_error:e}"
        );
        let final_totals = flow.totals();
        for i in 0..5 {
            assert!((final_totals.species_masses[i] - initial.species_masses[i]).abs() < 2e-13);
        }
        assert!((final_totals.energy - initial.energy).abs() < 2e-12);
    }
}
