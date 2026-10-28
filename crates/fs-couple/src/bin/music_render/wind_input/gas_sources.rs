//! Gas packs and solid/radiation artifacts share one admitted source snapshot.
use super::{attachments::AttachedSources, options, run_with_prepared_sources};
use super::super::thermochemical;
use std::path::Path;

pub(super) fn run_with_gas_options(
    args: &[String],
    gas_options: &thermochemical::Options,
    required_schema: Option<&'static str>,
) -> Result<(), String> {
    let (input, _, _, _, _, _, self_contained, attached_names) = options(args)?;
    if gas_options.record.is_some() {
        return Err(if self_contained {
            "--self-contained refuses external gas source files; wind gas records belong in the performance file"
        } else {
            "wind gas records belong in the performance file"
        }.into());
    }
    let attached = if self_contained {
        let parent = Path::new(input).parent().unwrap_or(Path::new("."));
        Some(AttachedSources::load(parent, &attached_names)?)
    } else { None };
    let loaded = if let Some(snapshot) = attached.as_ref() {
        gas_options.load_with(|path| {
            let name = path.to_str().ok_or_else(||
                "attached gas source paths must be UTF-8".to_string())?;
            snapshot.get(name).map(<[u8]>::to_vec).ok_or_else(|| format!(
                "--self-contained refuses external gas source files: {name:?} was not attached"
            ))
        })?
    } else {
        gas_options.load()?
    };
    let bindings = loaded.bindings();
    let registry = loaded.registry(&bindings);
    run_with_prepared_sources(args, registry.as_ref(), attached, required_schema)
}
