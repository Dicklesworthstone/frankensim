//! Source-bound film declarations for the existing material shell compiler.
//! Thickness, tensor-lift assumptions and spectral approximation limits are
//! authored inputs; the published source supplies the density and spectrum.
use super::*;
use fs_couple::shear_shell::{CoatingSide, ResolvedShearShell, ShellCoating};
use fs_material::coating::ResolvedCoating;
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

pub(super) const MAX_COATINGS: usize = 64;

#[derive(Clone, Debug)]
pub(super) struct FilmInput {
    pack: String,
    side: CoatingSide,
    thickness_m: f64,
    reference_hz: f64,
    poisson_ratio: f64,
    strain_limit: f64,
    relative_humidity: f64,
    storage_tolerance: f64,
}

impl FilmInput {
    /// PACK SIDE THICKNESS_M REFERENCE_HZ NU STRAIN_LIMIT RH STORAGE_TOL.
    pub(super) fn parse(values: &[&str]) -> Result<Self> {
        if values.len() != 8 {
            return Err(
                "shell coating needs eight material, geometry and approximation fields".into(),
            );
        }
        let pack = values[0];
        if pack.is_empty()
            || pack == "@body"
            || pack.chars().any(char::is_control)
            || Path::new(pack).is_absolute()
            || pack.contains('\\')
            || Path::new(pack)
                .components()
                .any(|c| !matches!(c, std::path::Component::Normal(_)))
        {
            return Err(
                "shell coating pack must be an explicit relative path without traversal".into(),
            );
        }
        let side = match values[1] {
            "positive" => CoatingSide::Positive,
            "negative" => CoatingSide::Negative,
            _ => return Err("shell coating side must be positive or negative".into()),
        };
        let film = Self {
            pack: pack.into(),
            side,
            thickness_m: scalar(values[2])?,
            reference_hz: scalar(values[3])?,
            poisson_ratio: scalar(values[4])?,
            strain_limit: scalar(values[5])?,
            relative_humidity: scalar(values[6])?,
            storage_tolerance: scalar(values[7])?,
        };
        if film.thickness_m <= 0.0
            || film.reference_hz <= 0.0
            || !(core::f64::consts::TAU * film.reference_hz).is_finite()
            || film.poisson_ratio <= -1.0
            || film.poisson_ratio >= 0.5
            || film.strain_limit <= 0.0
            || !(0.0..=1.0).contains(&film.relative_humidity)
            || !(0.0..1.0).contains(&film.storage_tolerance)
        {
            return Err("shell coating requires positive thickness/frequency/strain, -1 < nu < 0.5, 0 <= RH <= 1 and storage tolerance in [0,1)".into());
        }
        Ok(film)
    }

    pub(super) fn write_row(&self, region: &str, out: &mut String) -> Result<()> {
        use std::fmt::Write as _;
        let side = match self.side {
            CoatingSide::Positive => "positive",
            CoatingSide::Negative => "negative",
        };
        writeln!(
            out,
            "coating {region} {} {side} {} {} {} {} {} {}",
            self.pack,
            self.thickness_m,
            self.reference_hz,
            self.poisson_ratio,
            self.strain_limit,
            self.relative_humidity,
            self.storage_tolerance
        )?;
        Ok(())
    }

    fn resolve(
        &self,
        region: &str,
        temperature: f64,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
    ) -> Result<ShellCoating> {
        load_pack(&self.pack, catalog, loader)?;
        let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        // These axis names and flags are part of the measured-film schema.
        // In particular the group density is NOT a sample-matched DMA density.
        let query = QueryPoint::new()
            .with_quantity(
                "temperature",
                semantic(QuantityKind::AbsoluteTemperature),
                temperature,
            )?
            .with_quantity(
                "frequency",
                semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
                self.reference_hz,
            )?
            .with_quantity(
                "conditioning_relative_humidity",
                QuantitySpec::dimensional(Dims::NONE),
                self.relative_humidity,
            )?
            .with_quantity(
                "source_fitted_response",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )?
            .with_quantity(
                "source_density_sample_matched",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )?;
        let material = resolve_material_state_point(
            &catalog[&self.pack].card,
            &query,
            &ResolvedCoating::requirements()?,
            Pick::SingleClaimOnly,
        )?;
        Ok(ShellCoating {
            region: region.into(),
            material: ResolvedCoating::from_material(
                material,
                self.poisson_ratio,
                self.strain_limit,
            )?,
            side: self.side,
            thickness_m: self.thickness_m,
            reference_omega_rad_s: core::f64::consts::TAU * self.reference_hz,
            max_storage_relative_drift: self.storage_tolerance,
        })
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct Inputs(Vec<(String, FilmInput)>);

impl Inputs {
    pub(super) fn parse(rows: &mut Records<'_>, regions: &[RegionInput]) -> Result<Self> {
        if rows
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next())
            != Some("coatings")
        {
            return Ok(Self::default());
        }
        let n = count(rows.row("coatings", 1)?[0], 0, MAX_COATINGS)?;
        let mut declarations = Vec::with_capacity(n);
        for _ in 0..n {
            let v = rows.row("coating", 9)?;
            if !regions.iter().any(|region| region.name == v[0]) {
                return Err("shell coating refers to an unknown material region".into());
            }
            declarations.push((v[0].into(), FilmInput::parse(&v[1..])?));
        }
        Ok(Self(declarations))
    }

    pub(super) fn bind(
        &self,
        specimen: ResolvedShearShell,
        spec: &input::Spec,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        gate: &CancelGate,
    ) -> Result<ResolvedShearShell> {
        if self.0.is_empty() {
            return Ok(specimen);
        }
        let mut coatings = Vec::with_capacity(self.0.len());
        for (region, film) in &self.0 {
            poll(gate)?;
            coatings.push(film.resolve(region, spec.temperature, catalog, loader)?);
        }
        Ok(specimen.with_coatings(coatings, gate)?)
    }
}
