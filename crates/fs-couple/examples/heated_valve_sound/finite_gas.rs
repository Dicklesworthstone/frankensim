//! Finite gas heating of the same mesh used by the valve mechanics.
//! Internal links are a declared coarse, two-point thermal network, not a
//! consistent tetrahedral heat-equation discretization on a nonorthogonal mesh.
use super::*;
use std::collections::BTreeMap;

use fs_couple::solid_enthalpy::gas::{GasConductionGeometry, SolidGasContact, SolidGasSolve};
use fs_material::gas_species::input::read_isochoric_parcel;
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
    resolve_material_state_point,
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::{ThermalLink, series_conductance_w_k};

pub(super) struct Options {
    pub gas: IsochoricGasParcel,
    pub gas_distance_m: f64,
    pub steps: usize,
}
impl Options {
    pub fn extract(
        args: &[String],
        sources: &thermochemical::Options,
    ) -> Result<(Vec<String>, Option<Self>)> {
        let mut remaining = Vec::new();
        let mut specification = None;
        let mut words = args.iter();
        while let Some(word) = words.next() {
            if word != "--finite-gas" {
                remaining.push(word.clone());
                continue;
            }
            if specification.is_some() {
                return Err("--finite-gas may be supplied only once".into());
            }
            let mut field = || {
                words
                    .next()
                    .ok_or("--finite-gas needs volume, gas-layer distance and steps")
            };
            let volume: f64 = field()?.parse()?;
            let distance: f64 = field()?.parse()?;
            let steps: usize = field()?.parse()?;
            if !volume.is_finite()
                || volume <= 0.0
                || volume > 1.0
                || !distance.is_finite()
                || distance <= 0.0
                || distance > 1.0
                || !(1..=1024).contains(&steps)
            {
                return Err("finite gas requires volume in (0,1] m3, layer distance in (0,1] m and 1..=1024 steps".into());
            }
            specification = Some((volume, distance, steps));
        }
        let Some((volume, gas_distance_m, steps)) = specification else {
            return Ok((remaining, None));
        };
        let loaded = sources.load()?;
        let bindings = loaded.bindings();
        let registry = loaded
            .registry(&bindings)
            .ok_or("finite gas requires explicit component, water and frozen-phase declarations")?;
        let path = sources
            .record
            .as_ref()
            .ok_or("--finite-gas requires --gas-record")?;
        let mut bytes = Vec::new();
        std::fs::File::open(path)?
            .take((1 << 20) + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() > 1 << 20 {
            return Err("gas record exceeds 1 MiB".into());
        }
        let mut words = std::str::from_utf8(&bytes)?.split_ascii_whitespace();
        let gas = read_isochoric_parcel(&mut words, &registry, volume)?;
        if words.next().is_some() {
            return Err("gas record contains trailing fields".into());
        }
        if !(280.0..=400.0).contains(&gas.state().state().temperature) {
            return Err("example gas inversion interval is 280..400 K".into());
        }
        Ok((
            remaining,
            Some(Self {
                gas,
                gas_distance_m,
                steps,
            }),
        ))
    }
}

pub(super) struct Transport {
    pub links: Vec<ThermalLink>,
    pub surface_areas_m2: Vec<f64>,
    contacts: Vec<SolidGasContact>,
}

fn centroid(mesh: &Mesh, element: usize) -> [f64; 3] {
    core::array::from_fn(|axis| {
        mesh.tets[element]
            .iter()
            .map(|&node| mesh.nodes[node][axis])
            .sum::<f64>()
            / 4.0
    })
}

/// A constant source claim is necessary because the existing internal thermal
/// links are fixed in this bounded example. Endpoint sampling cannot prove that
/// an arbitrary temperature curve is constant, so non-scalar claims refuse.
fn constant_conductivity(card: &MaterialCard) -> Result<f64> {
    let requirement = ScalarPropertyRequirement::try_new(
        "thermal_conductivity",
        THERMAL_CONDUCTIVITY_DIMS,
        ScalarAdmissibility::StrictlyPositive,
    )?;
    let lower = resolve_material_state_point(
        card,
        &point(280.0)?,
        &[requirement.clone()],
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let upper = resolve_material_state_point(
        card,
        &point(400.0)?,
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let a = lower
        .property_by_key(&key)
        .ok_or("missing conductivity snapshot")?
        .answer();
    let b = upper
        .property_by_key(&key)
        .ok_or("missing conductivity snapshot")?
        .answer();
    if a.receipt.selected != b.receipt.selected {
        return Err("one constant conductivity source must span the entire solid interval".into());
    }
    let claim = card
        .claims()
        .claim(a.receipt.selected)
        .ok_or("selected conductivity claim missing")?;
    match (&claim.value, claim.interpolation) {
        (PropertyValue::Scalar { value, .. }, InterpolationPolicy::ConstantWithinValidity) => Ok(*value),
        _ => Err("finite-gas valve internal links require an explicit constant scalar conductivity claim".into()),
    }
}

pub(super) fn transport(
    mesh: &Mesh,
    geometry: &TetLoadGeometry3<'_>,
    card: &MaterialCard,
    options: &Options,
    cx: &Cx<'_>,
) -> Result<Transport> {
    let k = constant_conductivity(card)?;
    let mut owners = BTreeMap::<[usize; 3], Vec<usize>>::new();
    for (element, tet) in mesh.tets.iter().enumerate() {
        for local in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            let mut key = local.map(|i| tet[i]);
            key.sort_unstable();
            owners.entry(key).or_default().push(element);
        }
    }
    let mut links = Vec::new();
    let mut surface_areas_m2 = vec![0.0; mesh.tets.len()];
    for (face, elements) in &owners {
        cx.checkpoint().map_err(|_| "heated valve cancelled")?;
        let [a, b, c] = face.map(|node| mesh.nodes[node]);
        let u: [f64; 3] = core::array::from_fn(|i| b[i] - a[i]);
        let v: [f64; 3] = core::array::from_fn(|i| c[i] - a[i]);
        let normal = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        let area = 0.5 * length;
        for &element in elements {
            surface_areas_m2[element] += area;
        }
        if let [first, second] = elements.as_slice() {
            let distance = |element| {
                let center = centroid(mesh, element);
                (0..3)
                    .map(|i| (center[i] - a[i]) * normal[i] / length)
                    .sum::<f64>()
                    .abs()
            };
            links.push(ThermalLink {
                id: links.len() as u64,
                a: CellId(*first as u64),
                b: CellId(*second as u64),
                area_m2: area,
                conductance_w_k: series_conductance_w_k(
                    area,
                    [distance(*first), distance(*second)],
                    [k, k],
                )?,
            });
        }
    }
    // Reserve only a floating-point summation allowance, preventing cancellation
    // of the four face areas from producing a negative nominal exposed area.
    // This extra <=16-epsilon area is insulated and never becomes a gas patch.
    for area in &mut surface_areas_m2 {
        *area *= 1.0 + 16.0 * f64::EPSILON;
    }
    let mut contacts = Vec::new();
    for (id, face) in geometry.faces().iter().enumerate() {
        // The clamped root is insulated. All other exterior faces see the same
        // well-mixed parcel; interior faces only exchange with their neighbor.
        if face.nodes.iter().all(|&node| mesh.nodes[node][0] == 0.0) {
            continue;
        }
        let element = owners[&face.nodes][0];
        let center = centroid(mesh, element);
        let distance = (0..3)
            .map(|i| (center[i] - face.centroid_m[i]) * face.outward_normal[i])
            .sum::<f64>()
            .abs();
        contacts.push(SolidGasContact {
            id: id as u64,
            cell: CellId(element as u64),
            geometry: GasConductionGeometry {
                area_m2: face.area_m2,
                gas_distance_m: options.gas_distance_m,
                cell_distance_m: distance,
            },
            conductivity: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
    }
    Ok(Transport {
        links,
        surface_areas_m2,
        contacts,
    })
}

pub(super) struct History {
    // t, Tmin, Tmax, Tgas, pgas, external Q, solid dH, gas dU, solid->gas Q, residual
    pub rows: Vec<[f64; 10]>,
    pub description: String,
}
impl History {
    pub fn totals(&self) -> [f64; 3] {
        let last = self.rows.last().expect("initial history row");
        [last[5], last[6] + last[7], last[9]]
    }
    pub fn write(&self, path: &str) -> Result<()> {
        let mut out = std::io::BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)?,
        );
        writeln!(out, "# {}", self.description)?;
        writeln!(
            out,
            "time_s,solid_min_K,solid_max_K,gas_K,gas_Pa,external_heat_J,solid_storage_J,gas_storage_J,solid_to_gas_J,combined_residual_J"
        )?;
        for row in &self.rows {
            for (index, value) in row.iter().enumerate() {
                if index != 0 {
                    write!(out, ",")?;
                }
                write!(out, "{value:.17e}")?;
            }
            writeln!(out)?;
        }
        out.flush()?;
        eprintln!("{}", self.description);
        Ok(())
    }
}

pub(super) struct Evolution<'t> {
    pub thermal: DepositionSimulation<'t>,
    pub gas: IsochoricGasParcel,
    pub history: History,
}

pub(super) fn advance<'t>(
    source: &SolidEnthalpy<'_>,
    mut thermal: DepositionSimulation<'t>,
    duration_s: f64,
    options: &Options,
    transport: &Transport,
    cx: &Cx<'_>,
) -> Result<Evolution<'t>> {
    let mut gas = options.gas.clone();
    let initial_solid: BTreeMap<_, _> = thermal
        .cells()
        .iter()
        .map(|(&id, cell)| (id, cell.specific_h_j_kg))
        .collect();
    let initial_gas = gas.state().specific_internal_energy_j_kg();
    let controls = SolidGasSolve {
        maximum_cells: 48,
        maximum_contacts: 48,
        maximum_trials: 4096,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-7,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-7,
        gas_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-7,
        maximum_gas_temperature_increment_k: 20.0,
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 400.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 49,
            ..NewtonKrylovConfig::default()
        },
    };
    let provenance =
        thermochemical::provenance(
            &[fs_material::gas_species::input::GasRecord::Thermochemical(
                gas.state().clone(),
            )],
        );
    let mut history = History {
        rows: vec![[
            0.0,
            300.0,
            300.0,
            gas.state().state().temperature,
            gas.state().state().pressure,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]],
        description: format!(
            "finite fixed-volume gas; volume_m3={}; gas_mass_kg={}; gas_distance_m={}; steps={}; solid_interval_K=280..400; thermal network: fixed scalar source conductivity, centroid-normal two-point links, insulated clamped root; heating and sound are sequential frozen-state experiments; gas_source={provenance}",
            gas.volume_m3(),
            gas.mass_kg(),
            options.gas_distance_m,
            options.steps
        ),
    };
    let mut external = 0.0;
    let mut exchange = 0.0;
    for step in 1..=options.steps {
        let target = duration_s * step as f64 / options.steps as f64;
        let candidate = source.step_with_isochoric_gas(
            &thermal,
            &gas,
            target,
            0.0,
            &transport.contacts,
            controls,
            &vec![0.0; thermal.cells().len() + 1],
            cx,
        )?;
        let (next_thermal, next_gas, _, report) = candidate.into_parts();
        external += report.external_thermal_heat_j;
        exchange += report.solid_to_gas_heat_j;
        let solid = next_thermal
            .cells()
            .iter()
            .map(|(&id, cell)| cell.mass_kg * (cell.specific_h_j_kg - initial_solid[&id]))
            .sum::<f64>();
        let gas_storage =
            next_gas.mass_kg() * (next_gas.state().specific_internal_energy_j_kg() - initial_gas);
        let residual = solid + gas_storage - external;
        if !residual.is_finite() || residual.abs() > step as f64 * controls.total_energy_tolerance_j
        {
            return Err("cumulative endpoint solid/gas energy does not close".into());
        }
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        for &id in next_thermal.cells().keys() {
            let t = next_thermal.temperature_k(id)?;
            minimum = minimum.min(t);
            maximum = maximum.max(t);
        }
        history.rows.push([
            target,
            minimum,
            maximum,
            next_gas.state().state().temperature,
            next_gas.state().state().pressure,
            external,
            solid,
            gas_storage,
            exchange,
            residual,
        ]);
        thermal = next_thermal;
        gas = next_gas;
    }
    Ok(Evolution {
        thermal,
        gas,
        history,
    })
}

#[cfg(test)]
#[allow(dead_code)]
#[path = "../../tests/material_exchange/gas.rs"]
mod gas_fixture;

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::thermochemical_gas::ThermochemicalGasModel;
    use fs_qty::QuantitySpec;
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};

    fn options(volume: f64) -> Options {
        let (species, transport) = gas_fixture::component("heated-valve-fixture", 0.028, 0.0);
        let model = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            )
            .unwrap(),
            vec![transport],
            None,
        )
        .unwrap();
        let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let point = QueryPoint::new()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                340.0,
            )
            .unwrap()
            .with_quantity("pressure", typed(QuantityKind::Pressure), 100000.0)
            .unwrap();
        Options {
            gas: IsochoricGasParcel::new(model, point, volume).unwrap(),
            gas_distance_m: 0.001,
            steps: 8,
        }
    }

    #[test]
    fn finite_reservoir_size_changes_spatial_heating_valve_and_actual_sound() {
        let card = demo(1.0).unwrap();
        let small = options(0.0001);
        let large = options(0.0004);
        let a = synthesize_with_gas(&card, 0.0, 10.0, Some(&small), &CancelGate::new()).unwrap();
        let b = synthesize_with_gas(&card, 0.0, 10.0, Some(&large), &CancelGate::new()).unwrap();
        for result in [&a, &b] {
            let history = result.gas_history.as_ref().unwrap();
            assert_eq!(history.rows.len(), 9);
            let last = history.rows.last().unwrap();
            assert!(last[1] > 300.0 && last[2] > last[1] + 1e-5);
            assert!(last[3] < 340.0 && last[3] > last[2]);
            assert_eq!(last[5], 0.0);
            assert!(last[6] > 0.0 && last[7] < 0.0 && last[8] < 0.0);
            assert!((last[6] + last[7]).abs() < 2e-6);
            assert!((last[6] + last[8]).abs() < 2e-6);
            assert!((last[4] / 100000.0 - last[3] / 340.0).abs() < 1e-12);
            assert_eq!(result.cold_mass, result.hot_mass);
            assert!(result.displacement_m > 0.0);
            assert!(result.hot.iter().any(|p| *p != 0.0));
            assert_ne!(result.cold, result.hot);
        }
        assert_eq!(a.cold, b.cold);
        assert_ne!(a.hot, b.hot);
        assert_ne!(a.hot_frequency, b.hot_frequency);
        assert!(
            b.gas_history.as_ref().unwrap().rows.last().unwrap()[6]
                > a.gas_history.as_ref().unwrap().rows.last().unwrap()[6]
        );
        let replay =
            synthesize_with_gas(&card, 0.0, 10.0, Some(&small), &CancelGate::new()).unwrap();
        assert_eq!(a.hot, replay.hot);
        assert_eq!(
            a.gas_history.unwrap().rows,
            replay.gas_history.unwrap().rows
        );
    }

    #[test]
    fn invalid_native_gas_input_and_cancelled_trajectory_return_no_comparison() {
        let sources = thermochemical::Options::default();
        for words in [
            vec!["--finite-gas", "0", ".001", "8"],
            vec!["--finite-gas", ".0001", "0", "8"],
            vec!["--finite-gas", ".0001", ".001", "0"],
            vec!["--finite-gas", ".0001", ".001", "8"],
        ] {
            let args = words.into_iter().map(String::from).collect::<Vec<_>>();
            assert!(Options::extract(&args, &sources).is_err());
        }
        let gate = CancelGate::new();
        gate.request();
        assert!(
            synthesize_with_gas(
                &demo(1.0).unwrap(),
                0.0,
                10.0,
                Some(&options(0.0001)),
                &gate
            )
            .is_err()
        );
    }
}
