//! Finite contact-heat profile for the original material-string renderer.
//! Cells refer to material packs already used by this assembly. Density, Cp,
//! conductivity and film constitutive properties are queried from those SAME
//! cards. The profile supplies physical geometry and numerical resolution only.
//! Outer boundaries are explicitly adiabatic; links use original series-Fourier
//! conduction, with conductivity re-queried at every nonlinear endpoint.

use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    LayerThermalCell, LayerThermalSolve, SeriesDriveHeatDestination,
};
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve,
    SolidLiquidPhase,
};
use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent, backward_euler::StepCoefficients,
    series_conductance_w_k,
};

pub const MAX_BYTES: usize = 65_536;
const MAX_CELLS: usize = 256;
const MAX_LINKS: usize = 1024;

#[derive(Clone, Debug)]
struct Cell {
    id: CellId,
    pack: String,
    temperatures: [f64; 3], // lower source bound, initial, upper source bound
    geometry: CellGeometry,
}
#[derive(Clone, Debug)]
struct Link {
    id: u64,
    a: CellId,
    b: CellId,
    area: f64,
    distances: [f64; 2],
}
#[derive(Clone, Debug)]
pub struct Profile {
    cells: Vec<Cell>,
    films: Vec<(usize, CellId)>,
    series: Vec<SeriesDriveHeatDestination>,
    links: Vec<Link>,
    controls: LayerThermalSolve,
    caloric_error_j_kg: f64,
}
struct Row<'a> {
    words: std::str::SplitAsciiWhitespace<'a>,
}
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str> {
        self.words
            .next()
            .ok_or_else(|| "incomplete contact-thermal record".into())
    }
    fn parse<T: std::str::FromStr>(&mut self) -> Result<T> {
        self.word()?
            .parse()
            .map_err(|_| "invalid contact-thermal value".into())
    }
    fn scalar(&mut self) -> Result<f64> {
        let x: f64 = self.parse()?;
        if !x.is_finite() {
            return Err("contact-thermal values must be finite".into());
        }
        Ok(x)
    }
    fn positive(&mut self) -> Result<f64> {
        let x = self.scalar()?;
        if x <= 0.0 {
            return Err("contact-thermal value must be positive".into());
        }
        Ok(x)
    }
    fn end(&mut self) -> Result<()> {
        if self.words.next().is_some() {
            return Err("extra contact-thermal field".into());
        }
        Ok(())
    }
}
struct Rows<'a> {
    lines: std::str::Lines<'a>,
}
impl<'a> Rows<'a> {
    fn row(&mut self, name: &str) -> Result<Row<'a>> {
        let mut row = Row {
            words: self
                .lines
                .next()
                .ok_or("missing contact-thermal record")?
                .split_ascii_whitespace(),
        };
        if row.word()? != name {
            return Err(format!("expected contact-thermal {name} record").into());
        }
        Ok(row)
    }
    fn count(&mut self, name: &str, max: usize) -> Result<usize> {
        let mut row = self.row(name)?;
        let n = row.parse()?;
        row.end()?;
        if n > max {
            return Err(format!("contact-thermal {name} exceeds {max}").into());
        }
        Ok(n)
    }
}
impl Profile {
    pub fn parse(bytes: &[u8]) -> Result<Self> {
        if bytes.len() > MAX_BYTES {
            return Err("contact-thermal profile exceeds 64 KiB".into());
        }
        let mut rows = Rows {
            lines: std::str::from_utf8(bytes)?.lines(),
        };
        rows.row("frankensim-contact-thermal-v1")?.end()?;
        let mut boundary = rows.row("boundary")?;
        if boundary.word()? != "adiabatic" {
            return Err("contact-thermal v1 requires an explicit adiabatic outer boundary".into());
        }
        boundary.end()?;
        let mut row = rows.row("solver")?;
        let trials: usize = row.parse()?;
        let iterations: usize = row.parse()?;
        let difference_step_j = row.positive()?;
        let cell_energy_tolerance_j = row.positive()?;
        let thermal_energy_tolerance_j = row.positive()?;
        let total_energy_tolerance_j = row.positive()?;
        let absolute_tolerance = row.positive()?;
        let relative_tolerance = row.positive()?;
        let linear_restart: usize = row.parse()?;
        let caloric_error_j_kg = row.positive()?;
        row.end()?;
        if !(1..=16_384).contains(&trials)
            || !(1..=256).contains(&iterations)
            || !(1..=MAX_CELLS).contains(&linear_restart)
        {
            return Err(
                "contact-thermal solver exceeds its bounded trial/iteration/Krylov domain".into(),
            );
        }
        let count = rows.count("cells", MAX_CELLS)?;
        if count == 0 {
            return Err("contact-thermal network needs finite cells".into());
        }
        let mut cells = Vec::<Cell>::with_capacity(count);
        for _ in 0..count {
            let mut row = rows.row("cell")?;
            let id = CellId(row.parse()?);
            let pack = row.word()?.to_owned();
            if pack.len() > 512
                || !Path::new(&pack)
                    .components()
                    .all(|c| matches!(c, std::path::Component::Normal(_)))
            {
                return Err(
                    "contact-thermal cell needs an original relative material-pack path".into(),
                );
            }
            let temperatures = [row.positive()?, row.positive()?, row.positive()?];
            let start_m = [row.scalar()?, row.scalar()?, row.scalar()?];
            let end_m = [row.scalar()?, row.scalar()?, row.scalar()?];
            let geometry = CellGeometry {
                start_m,
                end_m,
                volume_m3: row.positive()?,
                surface_area_m2: row.positive()?,
            };
            row.end()?;
            if cells.iter().any(|c| c.id == id)
                || !(temperatures[0] < temperatures[1] && temperatures[1] < temperatures[2])
            {
                return Err("contact-thermal cells need distinct IDs and an initial temperature strictly inside source bounds".into());
            }
            cells.push(Cell {
                id,
                pack,
                temperatures,
                geometry,
            });
        }
        let mut films = Vec::new();
        for _ in 0..rows.count("films", 64)? {
            let mut row = rows.row("film")?;
            let port = row.parse()?;
            let cell = CellId(row.parse()?);
            row.end()?;
            if films.iter().any(|&(p, c)| p == port || c == cell)
                || !cells.iter().any(|c| c.id == cell)
            {
                return Err("each film needs one distinct physical port and existing cell".into());
            }
            films.push((port, cell));
        }
        let mut series = Vec::<SeriesDriveHeatDestination>::new();
        for _ in 0..rows.count("series", 64)? {
            let mut row = rows.row("series-heat")?;
            let route = SeriesDriveHeatDestination {
                port: row.parse()?,
                cell: CellId(row.parse()?),
            };
            row.end()?;
            if series.iter().any(|r| r.port == route.port)
                || !cells.iter().any(|c| c.id == route.cell)
            {
                return Err("series heat needs distinct physical drives and existing cells".into());
            }
            series.push(route);
        }
        if films.is_empty() && series.is_empty() {
            return Err("contact-thermal profile has no physical heat route".into());
        }
        let mut links = Vec::<Link>::new();
        for _ in 0..rows.count("links", MAX_LINKS)? {
            let mut row = rows.row("link")?;
            let link = Link {
                id: row.parse()?,
                a: CellId(row.parse()?),
                b: CellId(row.parse()?),
                area: row.positive()?,
                distances: [row.positive()?, row.positive()?],
            };
            row.end()?;
            if link.a == link.b
                || links.iter().any(|l| l.id == link.id)
                || !cells.iter().any(|c| c.id == link.a)
                || !cells.iter().any(|c| c.id == link.b)
            {
                return Err(
                    "contact-thermal links need distinct IDs and two different existing cells"
                        .into(),
                );
            }
            links.push(link);
        }
        if rows.lines.next().is_some() {
            return Err("trailing contact-thermal records".into());
        }
        let newton = NewtonKrylovConfig {
            absolute_tolerance,
            relative_tolerance,
            linear_restart,
            ..NewtonKrylovConfig::default()
        };
        newton.validate()?;
        Ok(Self {
            cells,
            films,
            series,
            links,
            caloric_error_j_kg,
            controls: LayerThermalSolve {
                maximum_cells: count,
                maximum_mechanical_trials: trials,
                maximum_newton_iterations: iterations,
                difference_step_j,
                cell_energy_tolerance_j,
                thermal_energy_tolerance_j,
                total_energy_tolerance_j,
                newton,
            },
        })
    }

    /// Compile thermal data from the SAME cached material objects as mechanics.
    /// A profile cannot silently load a different revision or substitute a
    /// constant for a missing Cp, density, conductivity or film property.
    pub fn prepare(&self, ready: &Ready, gate: &CancelGate) -> Result<Prepared> {
        poll(gate)?;
        let mut sources = Vec::with_capacity(self.cells.len());
        for cell in &self.cells {
            poll(gate)?;
            let card = ready.catalog.get(&cell.pack)
                .ok_or("contact-thermal cell must reference a material pack already used by this assembly")?.card.clone();
            let point = |t| material_inputs::query(&card, t, ready.spec.omega);
            let [low, initial, high] = cell.temperatures;
            let requirement = [
                Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?,
                Requirement::try_new(
                    "thermal_conductivity",
                    THERMAL_CONDUCTIVITY_DIMS,
                    Admit::NonNegative,
                )?,
            ];
            let state = resolve_material_state_point(
                &card,
                &point(initial)?,
                &requirement,
                Pick::SingleClaimOnly,
            )?;
            let cp = PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS);
            let density = PropertyKey::new("density", Density::DIMS);
            let span = |a, b| -> Result<ResolvedHeatCapacityCurve> {
                Ok(ResolvedHeatCapacityCurve::try_new(
                    &card,
                    HeatCapacitySourceQuery {
                        heat_capacity: &cp,
                        density: &density,
                        temperature_axis: "T",
                        lower: &point(a)?,
                        upper: &point(b)?,
                        selection: SelectionPolicy::SingleClaimOnly,
                    },
                    SolidLiquidPhase::Solid,
                    0.0,
                    self.caloric_error_j_kg,
                )?)
            };
            // Retain the initial temperature as an EXACT node with h=0.
            // Only the enthalpy datum shifts; native Cp integration is unchanged.
            let cold = span(low, initial)?;
            let hot = span(initial, high)?;
            let offset = cold
                .curve()
                .knots()
                .last()
                .ok_or("empty source caloric chart")?
                .specific_enthalpy_j_kg;
            let mut knots = cold.curve().knots().to_vec();
            for knot in &mut knots {
                knot.specific_enthalpy_j_kg -= offset;
            }
            knots.extend_from_slice(&hot.curve().knots()[1..]);
            let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
                card.content_hash(),
                SolidLiquidPhase::Solid,
                knots,
            )?;
            let rho = state
                .property("density")
                .ok_or("missing source density")?
                .value_si();
            if curve
                .knots()
                .iter()
                .any(|k| k.bulk_density_kg_m3.to_bits() != rho.to_bits())
            {
                return Err("fixed contact-thermal cells require one constant reference density over the full source interval".into());
            }
            sources.push(CellSource {
                card,
                state,
                curve,
                source_curves: [cold, hot],
                density: rho,
            });
        }
        let actual_films = (0..ready.spec.bow_contact_count())
            .filter(|&port| ready.performance.bank().bow_shear_layer(port).is_some())
            .count();
        if actual_films != self.films.len() {
            return Err("contact-thermal profile must bind every original material film".into());
        }
        for &(port, id) in &self.films {
            let index = self
                .cells
                .iter()
                .position(|c| c.id == id)
                .ok_or("missing film cell")?;
            let cell = &self.cells[index];
            let source = &sources[index];
            let layer = ready
                .performance
                .bank()
                .bow_shear_layer(port)
                .ok_or("film names an absent material shear port")?;
            if layer.resolved().card_identity() != source.card.content_hash()
                || layer.temperature_kelvin().to_bits() != cell.temperatures[1].to_bits()
                || layer.resolved().property("density").is_none()
                || (layer.area_m2() * layer.thickness_m() - cell.geometry.volume_m3).abs()
                    > 32.0 * f64::EPSILON * cell.geometry.volume_m3
            {
                return Err("film cell must retain its exact material, initial temperature, density and layer volume".into());
            }
            layer.thermal_stiffness_slope_n_m_k(
                &source.card,
                "T",
                (cell.temperatures[0], cell.temperatures[2]),
            )?;
        }
        poll(gate)?;
        let prepared = Prepared {
            profile: self.clone(),
            sources,
        };
        // Validate the source masses, geometry and initial conductances before
        // creating any output. This processes only birth/contact events at t=0.
        prepared.network(gate)?;
        Ok(prepared)
    }
}

struct CellSource {
    card: MaterialCard,
    state: ResolvedMaterialStatePoint,
    curve: EquilibriumEnthalpyPhaseCurve,
    source_curves: [ResolvedHeatCapacityCurve; 2],
    density: f64,
}
pub struct Prepared {
    profile: Profile,
    sources: Vec<CellSource>,
}
impl Prepared {
    fn conductance(&self, link: &Link, temperatures: [f64; 2]) -> Result<f64> {
        let mut k = [0.0; 2];
        for (i, id) in [link.a, link.b].iter().enumerate() {
            let index = self
                .profile
                .cells
                .iter()
                .position(|c| c.id == *id)
                .ok_or("missing thermal link cell")?;
            let source = &self.sources[index];
            let query = source.state.query().with("T", temperatures[i])?;
            let state = source.state.at_query_point(&source.card, &query)?;
            k[i] = state
                .property("thermal_conductivity")
                .ok_or("missing source conductivity")?
                .value_si();
        }
        Ok(series_conductance_w_k(link.area, link.distances, k)?)
    }
    fn network(&self, gate: &CancelGate) -> Result<DepositionSimulation<'_>> {
        let mut events = Vec::new();
        for (index, (cell, source)) in self.profile.cells.iter().zip(&self.sources).enumerate() {
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: cell.id,
                    segment_id: cell.id.0,
                    body_id: cell.id.0,
                    material_id: MaterialId(index as u64),
                    mass_kg: source.density * cell.geometry.volume_m3,
                    specific_h_j_kg: 0.0,
                    geometry: cell.geometry,
                    boundary: Boundary {
                        ambient_k: cell.temperatures[1],
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
        for link in &self.profile.links {
            let temperature = |id| {
                self.profile
                    .cells
                    .iter()
                    .find(|c| c.id == id)
                    .expect("admitted link")
                    .temperatures[1]
            };
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: link.id,
                    a: link.a,
                    b: link.b,
                    area_m2: link.area,
                    conductance_w_k: self
                        .conductance(link, [temperature(link.a), temperature(link.b)])?,
                }),
            });
        }
        let laws: Vec<_> = self
            .sources
            .iter()
            .enumerate()
            .map(|(i, source)| (MaterialId(i as u64), &source.curve as &dyn ThermalLaw))
            .collect();
        let network = DepositionSimulation::new(
            events,
            &laws,
            DepositionBudget {
                max_cells: self.profile.cells.len(),
                max_links: self.profile.links.len(),
                max_events: self.profile.cells.len() + self.profile.links.len(),
                max_substeps: 1,
            },
        )?
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: self.profile.controls.thermal_energy_tolerance_j,
            },
            || {
                if gate.is_requested() {
                    Err(fs_thermal::enthalpy::DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?
        .0;
        Ok(network)
    }
    fn resolve(&self, network: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>> {
        let mut conductances_w_k = BTreeMap::new();
        for link in &self.profile.links {
            conductances_w_k.insert(
                link.id,
                self.conductance(
                    link,
                    [
                        network.temperature_k(link.a)?,
                        network.temperature_k(link.b)?,
                    ],
                )?,
            );
        }
        Ok(StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k,
            payload: (),
        })
    }
    pub fn describe(&self) {
        eprintln!(
            "contact_thermal=source-bound finite solid cells; adiabatic exterior; fixed compressed films; original sample clock"
        );
        for (cell, source) in self.profile.cells.iter().zip(&self.sources) {
            eprintln!(
                "contact_thermal_cell={}; pack={:?}; card={}; caloric={}; source_caloric=[{},{}]; mass_kg={}; initial_T_K={}",
                cell.id.0,
                cell.pack,
                source.card.content_hash(),
                source.curve.identity(),
                source.source_curves[0].curve().identity(),
                source.source_curves[1].curve().identity(),
                source.density * cell.geometry.volume_m3,
                cell.temperatures[1]
            );
        }
    }
    /// Use the original scheduled pressure owner and PCM writer. Optional bow
    /// tracing observes accepted steps; it never advances a second renderer.
    pub fn write<W: Write + Seek>(
        &self,
        ready: &mut Ready,
        output: W,
        trace: Option<&mut dyn Write>,
        gate: &CancelGate,
    ) -> Result<(W, WavStreamSummary)> {
        self.write_observed(ready, output, trace, false, gate)
    }

    /// Thermal and mechanical state advance together at the source clock.
    /// Only their accepted pressure is decimated; traces retain every source step.
    pub fn write_observed<W: Write + Seek>(
        &self,
        ready: &mut Ready,
        output: W,
        trace: Option<&mut dyn Write>,
        decimate: bool,
        gate: &CancelGate,
    ) -> Result<(W, WavStreamSummary)> {
        let network = self.network(gate)?;
        let bindings = self
            .profile
            .films
            .iter()
            .map(|&(port, cell)| {
                let index = self
                    .profile
                    .cells
                    .iter()
                    .position(|c| c.id == cell)
                    .expect("admitted film");
                LayerThermalCell {
                    port,
                    cell,
                    material: &self.sources[index].card,
                    caloric: &self.sources[index].curve,
                    temperature_axis: "T",
                }
            })
            .collect();
        let mut performance = ready.performance.with_contact_thermal_network(
            network,
            bindings,
            self.profile.series.clone(),
            self.profile.controls,
            gate,
            |endpoint| self.resolve(endpoint),
        )?;
        let result = match trace {
            Some(trace) => write_traced_pressure(
                &mut performance, output, trace, &ready.spec, decimate, gate,
                |performance| performance.last_frame().map(|frame| &frame.mechanics),
            ),
            None => write_pressure(&mut performance, output, &ready.spec, decimate, gate),
        }?;
        for cell in performance.thermal().cells().values() {
            eprintln!(
                "contact_thermal_final_cell={}; T_K={}; stored_H_J={}",
                cell.id.0,
                performance.thermal().temperature_k(cell.id)?,
                cell.mass_kg * cell.specific_h_j_kg
            );
        }
        Ok(result)
    }
}
