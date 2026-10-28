//! Explicit transient SGS admission and output for the cooling scene reader.

use fs_lbm::conjugate::{FvSgsConfig, turbulence::TURBULENT_PRANDTL};
use fs_lbm::sgs3::{SgsModel, SgsUse};

use super::{
    FvFlow, J, Result, Scene, Transient, Turbulence, VoxelDomain, bad, num, number, quote,
};
use std::fmt::Write as _;

pub(super) const HELP: &str = "\nExploratory transient SGS: solver.sgs requires model (\"smagorinsky\" or\n\"wale\"), coefficient, filter_width_factor, maximum_eddy_to_molecular and\nusage: \"exploration\". It requires transient.flow: \"unsteady\" and marched\nenergy. Filter width is factor * cell_volume^(1/3); the ratio cap refuses\nrather than clips. Do not combine with solver.turbulence: \"lvel\",\nsteady-on-mean-flow energy or evolving fluid material properties. The same\nclosure drives momentum and heat (Pr_t=0.9). Results report the model,\nfilter range and final eddy-viscosity/conductivity diagnostics; no turbulent\nDNS or acceptance-render qualification is implied. Example:\n  frankensim --json cooling-cht examples/cooling-cht/sgs-heated-obstacle.json\n";

pub(super) fn parse(
    value: Option<&J>,
    turbulence: Turbulence,
    transient: Option<&Transient>,
    evolving_materials: bool,
) -> Result<Option<FvSgsConfig>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let object = value
        .as_object()
        .ok_or_else(|| bad("solver.sgs must be an object"))?;
    for (key, _) in object {
        if ![
            "model",
            "coefficient",
            "filter_width_factor",
            "maximum_eddy_to_molecular",
            "usage",
        ]
        .contains(&key.as_str())
        {
            return Err(bad(format!("unsupported solver.sgs field {key:?}")));
        }
    }
    let coefficient = number(value, "coefficient", "solver.sgs")?;
    let model = match value.str_field("model") {
        Some("smagorinsky") => SgsModel::Smagorinsky { cs: coefficient },
        Some("wale") => SgsModel::Wale { cw: coefficient },
        _ => return Err(bad("solver.sgs.model must be \"smagorinsky\" or \"wale\"")),
    };
    let usage = match value.str_field("usage") {
        Some("exploration") => SgsUse::Exploration,
        Some("acceptance-render") => SgsUse::AcceptanceRender,
        _ => {
            return Err(bad(
                "solver.sgs.usage must explicitly declare \"exploration\"",
            ));
        }
    };
    let config = FvSgsConfig {
        model,
        filter_width_factor: number(value, "filter_width_factor", "solver.sgs")?,
        maximum_eddy_to_molecular: number(value, "maximum_eddy_to_molecular", "solver.sgs")?,
        usage,
    };
    config
        .validate()
        .map_err(|e| bad(format!("solver.sgs: {e}")))?;
    if turbulence != Turbulence::Laminar {
        return Err(bad(
            "solver.sgs cannot be combined with solver.turbulence: lvel",
        ));
    }
    if evolving_materials {
        return Err(bad(
            "solver.sgs requires constant fluid properties; the evolving-material route is not an SGS owner",
        ));
    }
    let options = transient.and_then(|t| t.unsteady.as_ref()).ok_or_else(|| {
        bad("solver.sgs requires transient.flow: \"unsteady\"; steady and frozen-flow paths do not advance SGS")
    })?;
    if options.mean_flow_energy {
        return Err(bad(
            "solver.sgs requires transient.energy: \"march\"; steady-on-mean-flow omits the resolved velocity-temperature correlation",
        ));
    }
    Ok(Some(config))
}

fn model(config: FvSgsConfig) -> (&'static str, f64) {
    match config.model {
        SgsModel::Smagorinsky { cs } => ("smagorinsky", cs),
        SgsModel::Wale { cw } => ("wale", cw),
    }
}

/// Final projected-field diagnostics, not time averages or resolved spectra.
pub(super) struct Summary {
    config: FvSgsConfig,
    min_filter_m: f64,
    max_filter_m: f64,
    max_eddy: f64,
    max_ratio: f64,
    mean_ratio: f64,
    max_conductivity: f64,
}

impl Summary {
    pub(super) fn new(scene: &Scene, domain: &VoxelDomain, flow: &FvFlow) -> Option<Self> {
        let config = scene.sgs?;
        let (mut volume, mut weighted_eddy, mut max_eddy, mut max_heat_eddy) =
            (0.0, 0.0, 0.0_f64, 0.0_f64);
        let (mut min_filter_m, mut max_filter_m) = (f64::INFINITY, 0.0_f64);
        for c in 0..domain.cell_count() {
            if !domain.is_fluid(c) {
                continue;
            }
            let cell_volume = domain.volume(c);
            volume += cell_volume;
            weighted_eddy += cell_volume * flow.eddy_viscosity_m2_s[c];
            max_eddy = max_eddy.max(flow.eddy_viscosity_m2_s[c]);
            max_heat_eddy = max_heat_eddy.max(flow.heat_eddy_viscosity_m2_s[c]);
            let width = config.filter_width_factor * fs_math::det::pow(cell_volume, 1.0 / 3.0);
            min_filter_m = min_filter_m.min(width);
            max_filter_m = max_filter_m.max(width);
        }
        let nu = scene.fluid.kinematic_viscosity_m2_s;
        Some(Self {
            config,
            min_filter_m,
            max_filter_m,
            max_eddy,
            max_ratio: max_eddy / nu,
            mean_ratio: weighted_eddy / volume / nu,
            max_conductivity: scene.fluid.volumetric_heat_capacity() * max_heat_eddy
                / TURBULENT_PRANDTL,
        })
    }

    /// Append members inside the existing flow object.
    pub(super) fn json(&self, out: &mut String) -> Result<()> {
        let (name, coefficient) = model(self.config);
        let _ = write!(
            out,
            ",\"turbulence\":{},\"max_eddy_viscosity_ratio\":{},\"sgs\":{{\"model\":{},\"coefficient\":{},\"filter_width_factor\":{},\"maximum_eddy_to_molecular\":{},\"usage\":\"exploration\",\"thermal_prandtl\":{},\"min_filter_width_m\":{},\"max_filter_width_m\":{},\"final_max_eddy_viscosity_m2_s\":{},\"final_volume_mean_eddy_viscosity_ratio\":{},\"final_max_turbulent_conductivity_w_m_k\":{}}}",
            quote(name),
            num(self.max_ratio)?,
            quote(name),
            num(coefficient)?,
            num(self.config.filter_width_factor)?,
            num(self.config.maximum_eddy_to_molecular)?,
            num(TURBULENT_PRANDTL)?,
            num(self.min_filter_m)?,
            num(self.max_filter_m)?,
            num(self.max_eddy)?,
            num(self.mean_ratio)?,
            num(self.max_conductivity)?
        );
        Ok(())
    }

    pub(super) fn text(&self, out: &mut String) {
        let (name, coefficient) = model(self.config);
        let _ = writeln!(
            out,
            "turbulence={name}\nsgs_usage=exploration\nsgs_coefficient={coefficient}\nsgs_filter_width_factor={}\nsgs_maximum_eddy_to_molecular={}\nsgs_min_filter_width_m={}\nsgs_max_filter_width_m={}\nsgs_thermal_prandtl={TURBULENT_PRANDTL}\nfinal_max_eddy_viscosity_m2_s={}\nmax_eddy_viscosity_ratio={}\nfinal_volume_mean_eddy_viscosity_ratio={}\nfinal_max_turbulent_conductivity_w_m_k={}",
            self.config.filter_width_factor,
            self.config.maximum_eddy_to_molecular,
            self.min_filter_m,
            self.max_filter_m,
            self.max_eddy,
            self.max_ratio,
            self.mean_ratio,
            self.max_conductivity
        );
    }
}
