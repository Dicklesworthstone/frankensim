//! End-to-end (`f85xj.11.7`): a component power-budget table ON DISK becomes
//! an admitted [`PowerMap`] and a conduction solve.
//!
//! ## Why this test spans two crates
//!
//! `fs-conduction` deliberately takes no `fs-io` dependency: file parsing is a
//! quarantine concern, not a physics one, and the production graph stays free
//! of the importer. `fs-io` is therefore a DEV-dependency here — layer-exempt,
//! out of the production cone, and enough to prove the two halves compose.
//! The seam being tested is the one `f85xj.5.9` left open: `5.9` admits
//! already-typed [`ComponentPower`] rows, and nothing turned a file into them.
//!
//! ## What the ingestion boundary is actually protecting
//!
//! `fs-conduction` has NO unit guard. `COMPONENT_POWER_DIMS` is defined and
//! re-exported and is read nowhere in the workspace, so
//! [`ComponentPower::new`] takes a bare `f64` and a milliwatt figure is
//! accepted in silence — the audit then balances perfectly, at 1/1000 of the
//! real power, and the solve returns a plausible temperature field that is
//! wrong by three orders of magnitude. The `#power-unit:` directive at the
//! fs-io boundary is the only place that error can be caught, and
//! [`the_declared_unit_moves_the_solved_temperature_by_exactly_its_factor`]
//! pins the size of what is being prevented rather than describing it.
//!
//! ## Gauntlet tiers
//!
//! G0 for the exact delivered/declared identity and the refusal surface; G3
//! for the unit-rescaling metamorphic relation (identical digits under a
//! different declared unit must move the solved temperature rise by exactly
//! the unit factor).

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule, solve,
};
use fs_conduction::{
    ComponentPower, ConductionError, ConductionSolution, DEFAULT_POWER_TOLERANCE, PowerAudit,
    PowerMap, PowerUncertainty,
};
use fs_io::power_table::{PowerTable, PowerUnit, parse_power_table};
use fs_rep_mesh::TetComplex;
use support::with_cx;

/// The tracked power-budget export. Read from disk rather than
/// `include_str!`-ed: the DONE-WHEN is a table ON DISK, and a compiled-in
/// string would not exercise the caller's read at all.
const BUDGET_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/data/component-power-budget.csv"
);

/// The board the fixture describes: 120 mm x 80 mm x 2 mm, edge-cooled.
const EXTENT: [f64; 3] = [0.12, 0.08, 0.02];
/// One cell per 20 mm in plane, two through thickness, so every footprint
/// boundary below lands exactly on a grid line.
const CELLS: [usize; 3] = [6, 4, 2];
/// Effective in-plane conductivity of a metal-core board, W/(m K).
const CONDUCTIVITY: f64 = 20.0;
/// The cold rail temperature, K.
const RAIL_K: f64 = 300.0;
/// Grid coordinates are exact multiples of the cell size, so containment
/// needs only enough slack to absorb binary64 representation of the bounds.
const GEOMETRIC_EPS: f64 = 1e-9;

/// The name -> region seam, stated explicitly because the bead requires it to
/// be explicit.
///
/// `fs_io::power_table` produces NAMES and WATTS and resolves nothing to
/// geometry; a parser that resolved names would need a mesh. Binding a name to
/// a vertex set belongs to the layer that owns persistent entity identity.
/// This fixture-local catalog stands in for that layer so the seam is
/// exercised here instead of assumed, and so its failure mode is pinned: an
/// unknown component REFUSES. Skipping it would be far worse than a crash —
/// the model would quietly lose real dissipation and still balance, because
/// the declared total would have to be re-derived from whatever survived.
struct FootprintCatalog {
    /// `(component, x-range, y-range)`; each footprint spans the full board
    /// thickness, which is the lumped assumption a board-level budget makes.
    entries: &'static [(&'static str, [f64; 2], [f64; 2])],
}

impl FootprintCatalog {
    /// The Rev C placement. Every footprint clears the `x = 0.12` cold rail,
    /// so no component's power lands on a Dirichlet row where it would be
    /// absorbed by the boundary rather than solved for.
    const fn rev_c() -> Self {
        Self {
            entries: &[
                ("cpu-die", [0.00, 0.04], [0.00, 0.04]),
                ("dram-stack", [0.06, 0.10], [0.00, 0.02]),
                ("pmic", [0.06, 0.08], [0.06, 0.08]),
                ("phy", [0.00, 0.04], [0.06, 0.08]),
            ],
        }
    }

    /// Resolve one component name to the vertices it heats.
    ///
    /// Refuses through [`ConductionError::ScenarioRow`] rather than a new
    /// error variant: `ConductionError` is not `#[non_exhaustive]`, so adding
    /// a variant would be a breaking change at every match site, and the
    /// region/what/fix shape is exactly what an unresolvable declaration is.
    fn resolve(&self, mesh: &ConductionMesh, name: &str) -> Result<Vec<usize>, ConductionError> {
        let Some(&(_, x, y)) = self.entries.iter().find(|(entry, _, _)| *entry == name) else {
            return Err(ConductionError::ScenarioRow {
                region: name.to_owned(),
                what: format!(
                    "the power table declares component {name:?}, which has no footprint in the placement catalog"
                ),
                fix: "place the component in the scenario, or correct the component name in the power table; \
                      an unplaced component cannot be dropped, because its dissipation is real"
                    .to_owned(),
            });
        };
        let vertices: Vec<usize> = mesh
            .positions()
            .iter()
            .enumerate()
            .filter(|(_, p)| in_range(p[0], x) && in_range(p[1], y))
            .map(|(index, _)| index)
            .collect();
        if vertices.is_empty() {
            return Err(ConductionError::ScenarioRow {
                region: name.to_owned(),
                what: format!("footprint for {name:?} contains no mesh vertex"),
                fix: "refine the mesh or widen the footprint so the component has somewhere to deposit its power"
                    .to_owned(),
            });
        }
        Ok(vertices)
    }

    /// Whether a position lies in the named component's footprint.
    fn contains(&self, name: &str, position: [f64; 3]) -> bool {
        self.entries
            .iter()
            .find(|(entry, _, _)| *entry == name)
            .is_some_and(|(_, x, y)| in_range(position[0], *x) && in_range(position[1], *y))
    }
}

fn in_range(value: f64, range: [f64; 2]) -> bool {
    value >= range[0] - GEOMETRIC_EPS && value <= range[1] + GEOMETRIC_EPS
}

fn board_mesh() -> ConductionMesh {
    let (complex, positions) = box_grid(CELLS, EXTENT);
    let complex = TetComplex::from_tets(positions.len(), complex.tets);
    ConductionMesh::new(complex, positions).expect("board mesh admits")
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            tolerance: 1e-14,
            max_iterations: 60_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    }
}

/// Lower a parsed table onto a mesh: the whole ingestion path in one place.
///
/// The declared total is carried from the file VERBATIM into
/// [`PowerMap::new`]. That is deliberate — re-deriving it from the rows would
/// make the map balance by construction and destroy the only check that
/// catches a transcription error in the spreadsheet.
fn lower(
    table: &PowerTable,
    mesh: &ConductionMesh,
    catalog: &FootprintCatalog,
) -> Result<PowerMap, ConductionError> {
    let declared_total_w = table
        .declared_total_w()
        .expect("this fixture declares a total; a table without one cannot be audited");
    let mut components = Vec::with_capacity(table.rows().len());
    for row in table.rows() {
        let uncertainty = match (row.half_width_w(), row.confidence()) {
            (Some(half_width), Some(confidence)) => PowerUncertainty::HalfWidth {
                half_width,
                confidence,
            },
            _ => PowerUncertainty::Unstated,
        };
        components.push(ComponentPower::new(
            row.name(),
            row.watts(),
            uncertainty,
            catalog.resolve(mesh, row.name())?,
        )?);
    }
    PowerMap::new(components, declared_total_w)
}

/// Read, parse, lower, and solve. Returns the audit alongside the solution
/// because `ConductionSolution` carries no power audit — the audit exists only
/// where the source field was built, and a caller who wants it must hold it.
fn ingest_and_solve(text: &str, mesh: &ConductionMesh) -> (PowerAudit, ConductionSolution) {
    let table = parse_power_table(text).expect("fixture is well formed");
    let map = lower(&table, mesh, &FootprintCatalog::rev_c()).expect("map admits");
    let (source, audit) = map
        .volumetric_source(mesh, DEFAULT_POWER_TOLERANCE)
        .expect("projection admits");

    let boundary = ThermalBoundaryBuilder::new(mesh)
        .region(
            "cold-rail",
            |face| on_box_face(face.centroid[0], EXTENT[0]),
            ThermalBc::dirichlet(RAIL_K).expect("rail temperature"),
        )
        .expect("rail region")
        .adiabatic_remainder()
        .finish()
        .expect("complete partition");
    let material = ConductivityModel::isotropic_declared(CONDUCTIVITY).expect("board material");

    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config(),
        )
        .expect("solve")
    });
    (audit, solution)
}

fn budget_text() -> String {
    std::fs::read_to_string(BUDGET_PATH).expect("tracked power-budget fixture is readable")
}

#[test]
fn a_power_budget_export_on_disk_becomes_an_admitted_power_map_and_a_solve() {
    let text = budget_text();
    let table = parse_power_table(&text).expect("fixture is well formed");

    // The receipt is the quarantine provenance: what unit the file declared,
    // what factor was applied, and how many rows were admitted. A downstream
    // reviewer reads this, not the parser's source.
    let receipt = table.receipt();
    assert_eq!(receipt.declared_unit(), PowerUnit::Watts);
    assert!((receipt.watts_factor() - 1.0).abs() < f64::EPSILON);
    assert_eq!(receipt.rows(), 4);
    assert_eq!(table.declared_total_w(), Some(18.4));

    let mesh = board_mesh();
    let map = lower(&table, &mesh, &FootprintCatalog::rev_c()).expect("map admits");

    // PowerMap sorts by name, so downstream order is a property of the model
    // and not of the row order in whichever sheet the budget was exported
    // from. The fixture is deliberately not in alphabetical order.
    let names: Vec<&str> = map.components().iter().map(ComponentPower::name).collect();
    assert_eq!(names, ["cpu-die", "dram-stack", "phy", "pmic"]);

    let (audit, solution) = ingest_and_solve(&text, &mesh);

    // The identity that makes ingestion trustworthy: what the operator will
    // actually inject equals what the file declared. `delivered_total_w` is
    // recomputed from the ASSEMBLED nodal field (Σ f_a w_a), so this is not
    // per-component bookkeeping agreeing with itself.
    assert!(
        (audit.delivered_total_w() - 18.4).abs() < 1e-12,
        "delivered {} W against a declared 18.4 W",
        audit.delivered_total_w()
    );
    assert!((audit.component_total_w() - audit.declared_total_w()).abs() < 1e-12);

    // Per-row wattage survives the parse bit for bit; the conversion factor
    // for watts is exactly 1.0, so nothing may perturb the digits.
    let cpu = audit
        .rows()
        .iter()
        .find(|row| row.name() == "cpu-die")
        .expect("cpu-die is audited");
    assert!((cpu.declared_w() - 11.2).abs() < f64::EPSILON);
    assert!((cpu.delivered_w() - 11.2).abs() < 1e-12);

    // `phy` states no interval, so the aggregate band is None. A total half
    // width computed from the three rows that DID state one would be a
    // fabricated bound over a partially specified table — the weakest
    // ingredient limits the claim.
    assert_eq!(audit.total_half_width_w(), None);

    // The audit is the logged artifact the bead asks for. Assert its content
    // rather than only printing it, so a silently emptied table fails here.
    let rendered = audit.render_table();
    for component in ["cpu-die", "dram-stack", "pmic", "phy"] {
        assert!(
            rendered.contains(component),
            "audit table omits {component}:\n{rendered}"
        );
    }
    assert!(rendered.contains("declared 18.4 W"), "{rendered}");
    println!("{rendered}");

    // The solve agrees with the audit independently: the assembled energy
    // balance recomputes the source integral from the discrete operator.
    let energy = &solution.report.energy;
    assert!(
        (energy.source_w - audit.delivered_total_w()).abs() < 1e-9,
        "assembled source {} W against audited {} W",
        energy.source_w,
        audit.delivered_total_w()
    );
    assert!(
        energy.relative_closure() < 1e-9,
        "energy closure {} is not tight",
        energy.relative_closure()
    );
    // All 18.4 W leaves through the one cold rail, so the Dirichlet rows must
    // absorb it: sign included, because a sign slip here would still balance.
    assert!(
        (energy.dirichlet_in_w + 18.4).abs() < 1e-9,
        "rail carries {} W",
        energy.dirichlet_in_w
    );

    // The physics followed the table. The largest dissipator sits farthest
    // from the rail, so the hot spot must land inside its footprint — an
    // assertion that fails if the name->region seam mis-binds, which balance
    // checks alone cannot catch.
    let catalog = FootprintCatalog::rev_c();
    let (hottest, peak) = solution
        .temperature
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .expect("mesh has vertices");
    assert!(
        catalog.contains("cpu-die", mesh.positions()[hottest]),
        "hot spot {peak} K at {:?} is outside the cpu-die footprint",
        mesh.positions()[hottest]
    );
    assert!(*peak > RAIL_K, "heated board must exceed its rail");
}

#[test]
fn the_declared_unit_moves_the_solved_temperature_by_exactly_its_factor() {
    // G3, and the reason the directive is mandatory. IDENTICAL DIGITS under a
    // different declared unit: nothing downstream can tell these apart,
    // because fs-conduction never checks a dimension. Both solve, both audit
    // clean, and they differ by exactly 1000x.
    let watts_text = budget_text();
    let milliwatts_text = watts_text.replace("#power-unit: W", "#power-unit: mW");
    assert_eq!(
        watts_text.matches("#power-unit: W").count(),
        1,
        "the fixture's unit directive must be substitutable exactly once, \
         or this test silently compares a file with itself"
    );
    assert_ne!(watts_text, milliwatts_text);
    assert_eq!(
        watts_text.replace("#power-unit: W", ""),
        milliwatts_text.replace("#power-unit: mW", ""),
        "only the unit directive may differ"
    );

    let mesh = board_mesh();
    let (watt_audit, watt_solution) = ingest_and_solve(&watts_text, &mesh);
    let (milliwatt_audit, milliwatt_solution) = ingest_and_solve(&milliwatts_text, &mesh);

    assert!(
        (milliwatt_audit.delivered_total_w() * 1000.0 - watt_audit.delivered_total_w()).abs()
            < 1e-9,
        "mW delivered {} W against W delivered {} W",
        milliwatt_audit.delivered_total_w(),
        watt_audit.delivered_total_w()
    );

    // The consequence, in the quantity an engineer actually reads. The
    // problem is linear in the source, so the temperature RISE above the rail
    // scales exactly with the unit factor: a mis-declared table does not
    // produce an obviously broken field, it produces a 69 K hot spot reported
    // as 0.07 K.
    let watt_rise = watt_solution
        .temperature
        .iter()
        .map(|t| t - RAIL_K)
        .fold(0.0f64, f64::max);
    let milliwatt_rise = milliwatt_solution
        .temperature
        .iter()
        .map(|t| t - RAIL_K)
        .fold(0.0f64, f64::max);
    assert!(
        watt_rise > 1.0,
        "watt case must be a real rise: {watt_rise}"
    );
    assert!(
        (milliwatt_rise * 1000.0 - watt_rise).abs() < 1e-6 * watt_rise,
        "peak rise {milliwatt_rise} K (mW) against {watt_rise} K (W) is not a factor of exactly 1000"
    );
    for (index, (hot, cold)) in watt_solution
        .temperature
        .iter()
        .zip(&milliwatt_solution.temperature)
        .enumerate()
    {
        let scaled = (cold - RAIL_K) * 1000.0;
        assert!(
            (scaled - (hot - RAIL_K)).abs() < 1e-6 * watt_rise,
            "vertex {index}: {scaled} K against {} K",
            hot - RAIL_K
        );
    }
}

#[test]
fn a_declared_total_that_disagrees_with_its_rows_is_refused_by_the_power_map_not_the_parser() {
    // The adversarial case the bead names. The parser must NOT correct or
    // reject a disagreeing total: it carries it through verbatim so that
    // PowerMap's existing balance refusal is the one that fires. A second
    // audit implementation at the parse boundary would drift from this one.
    let text = budget_text().replace("#declared-total: 18.4", "#declared-total: 12.4");
    let table = parse_power_table(&text).expect("a disagreeing total is not a PARSE fault");
    assert_eq!(table.declared_total_w(), Some(12.4));
    assert!((table.row_total_w() - 18.4).abs() < 1e-12);

    let mesh = board_mesh();
    // The map itself admits: `PowerMap::new` checks only that the declared
    // total is finite and non-negative. The balance is a projection-time
    // property, because it is the projection that must deposit the power.
    let map = lower(&table, &mesh, &FootprintCatalog::rev_c()).expect("rows and total are legal");
    let refusal = map
        .volumetric_source(&mesh, DEFAULT_POWER_TOLERANCE)
        .expect_err("a map that does not add up must refuse");

    let ConductionError::ScenarioRow { region, what, .. } = &refusal else {
        panic!("expected a scenario-row refusal, got {refusal:?}");
    };
    assert_eq!(region, "<power-map>");
    assert!(what.contains("12.4"), "{what}");
    assert!(what.contains("18.4"), "{what}");
    // The refusal carries the whole per-component table, so the reviewer can
    // see which row is the transcription error without rerunning anything.
    assert!(what.contains("cpu-die"), "{what}");
}

#[test]
fn a_component_with_no_footprint_refuses_rather_than_dropping_its_power() {
    // The name->region seam's failure mode. This is the dangerous one: an
    // unresolved name is trivially "handled" by skipping the row, and the
    // result balances if the total is re-derived. Refusing is the only
    // behaviour that keeps the declared total meaningful.
    let text = budget_text().replace("pmic,2.1", "unplaced-heater,2.1");
    let table = parse_power_table(&text).expect("renaming a component is not a parse fault");

    let mesh = board_mesh();
    let refusal = lower(&table, &mesh, &FootprintCatalog::rev_c())
        .expect_err("an unplaced component must refuse");
    let ConductionError::ScenarioRow { region, what, .. } = &refusal else {
        panic!("expected a scenario-row refusal, got {refusal:?}");
    };
    assert_eq!(region, "unplaced-heater");
    assert!(what.contains("no footprint"), "{what}");

    // And the reason it matters: had the row been skipped, the remaining
    // components would sum to 16.3 W against a declared 18.4 W. That is a
    // 2.1 W hole — 11% of the budget — which the balance check would have
    // caught only because the total was NOT re-derived.
    let survivors: f64 = table
        .rows()
        .iter()
        .filter(|row| row.name() != "unplaced-heater")
        .map(fs_io::power_table::PowerTableRow::watts)
        .sum();
    assert!((survivors - 16.3).abs() < 1e-12, "{survivors}");
}
