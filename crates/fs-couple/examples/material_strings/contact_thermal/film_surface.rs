//! Material-attached film and sliding counterface on the original thermal owner.
use super::*;
use super::dry::{Strip, parse_surface};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::FrictionPortFrame;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::film_surface::RetainedFilmHeat;

#[derive(Clone, Debug)]
pub(super) struct Route {
    port: usize,
    direction: f64,
    film: Strip,
    surface: Vec<Strip>,
}

pub(super) fn parse(rows: &mut Rows<'_>, cells: &[Cell],
    films: &[(usize, CellId)]) -> Result<Vec<Route>> {
    let mut routes: Vec<Route> = Vec::new();
    for _ in 0..rows.count("film-surfaces", 64)? {
        let mut row = rows.row("film-surface")?;
        let port = row.parse()?;
        let direction = row.scalar()?;
        row.end()?;
        let cell = films.iter().find(|&&(p, _)| p == port)
            .map(|&(_, cell)| cell).ok_or("film surface requires its original retained film route")?;
        if routes.iter().any(|route| route.port == port) || ![-1.0, 1.0].contains(&direction) {
            return Err("film surfaces require distinct original ports and an explicit +1 or -1 tangent orientation".into());
        }
        let mut row = rows.row("film-footprint")?;
        let film = Strip { cell, interval: [row.scalar()?, row.scalar()?],
            width: row.positive()?, half_depth: row.positive()? };
        row.end()?;
        if film.interval[0] >= film.interval[1] {
            return Err("retained film footprint needs a positive material-coordinate length".into());
        }
        let mut occupied = std::collections::BTreeSet::from([cell]);
        let surface = parse_surface(rows, cells, &mut occupied)?;
        routes.push(Route { port, direction, film, surface });
    }
    Ok(routes)
}

pub(super) struct PreparedRoute<'a> {
    heat: RetainedFilmHeat<'a>,
    interval: [f64; 2],
    width: f64,
    direction: f64,
    bounds: [f64; 2],
}

impl PreparedRoute<'_> {
    fn footprint(&self, ports: &[FrictionPortFrame]) -> Result<[f64; 2]> {
        let layer = ports.get(self.heat.port()).and_then(|port| port.shear_layer.as_ref())
            .ok_or("retained film heat requires its actual candidate shear history")?;
        // The film remains on its original patch. Positive relative motion
        // brings smaller counterface coordinates under a positively oriented patch.
        let shift = -self.direction * (layer.elastic_slip_m + layer.plastic_slip_m);
        let mut interval = self.interval.map(|x| x + shift);
        if !shift.is_finite() || interval.iter().any(|x| !x.is_finite()) {
            return Err("counterface material displacement is unrepresentable".into());
        }
        let roundoff = 8.0 * f64::EPSILON * shift.abs()
            .max(self.interval[0].abs()).max(self.interval[1].abs());
        for edge in 0..2 {
            if (interval[edge] - self.bounds[edge]).abs() <= roundoff {
                let correction = self.bounds[edge] - interval[edge];
                interval = interval.map(|x| x + correction);
                break;
            }
        }
        Ok(interval)
    }
}

impl Prepared {
    pub(super) fn film_routes<'a>(&'a self, ready: &Ready,
        network: &DepositionSimulation<'_>) -> Result<Vec<PreparedRoute<'a>>> {
        let mut routes = Vec::with_capacity(self.profile.film_surfaces.len());
        for route in &self.profile.film_surfaces {
            let layer = ready.performance.bank().bow_shear_layer(route.port)
                .ok_or("counterface heat lost its original material film")?;
            let film = self.thermal_strip(&route.film)?;
            let binding = LayerThermalCell { port: route.port, cell: route.film.cell,
                material: film.material.material, caloric: film.caloric, temperature_axis: "T" };
            let surface = route.surface.iter().map(|strip| self.thermal_strip(strip))
                .collect::<Result<Vec<_>>>()?;
            let heat = RetainedFilmHeat::try_new(layer, binding, film, surface, network,
                self.profile.controls.maximum_cells)?;
            // Prove original coverage before output is opened. Every candidate
            // repeats this check at its actual displaced material coordinates.
            heat.heat_rates_w(network, route.film.interval, route.film.width)?;
            routes.push(PreparedRoute { heat, interval: route.film.interval,
                width: route.film.width, direction: route.direction,
                bounds: [route.surface[0].interval[0], route.surface.last().unwrap().interval[1]] });
        }
        Ok(routes)
    }

    pub(super) fn resolve_with_films(&self, endpoint: &DepositionSimulation<'_>,
        ports: &[FrictionPortFrame], routes: &[PreparedRoute<'_>]) -> Result<StepCoefficients<()>> {
        let mut coefficients = self.resolve(endpoint)?;
        let mut internal_residual_w = 0.0;
        for route in routes {
            let report = route.heat.heat_rates_w(endpoint, route.footprint(ports)?, route.width)?;
            internal_residual_w += report.internal_transfer_residual_w;
            for (cell, watts) in report.sources_w {
                let prior = coefficients.sources_w.entry(cell).or_insert(0.0);
                let sum = *prior + watts;
                let increment = sum - *prior;
                if !sum.is_finite() || !increment.is_finite()
                    || watts != 0.0 && (increment == 0.0 || increment.is_sign_positive() != watts.is_sign_positive())
                    || ((increment - watts) * self.sample_period_s).abs()
                        > self.profile.controls.total_energy_tolerance_j {
                    return Err("retained-film heat vanished in the original receiving-cell source sum".into());
                }
                internal_residual_w += increment - watts;
                *prior = sum;
            }
        }
        if !internal_residual_w.is_finite()
            || (self.sample_period_s * internal_residual_w).abs()
                > self.profile.controls.total_energy_tolerance_j {
            return Err("retained-film counterface transfer does not close its original energy account".into());
        }
        Ok(coefficients)
    }
}
