import { el } from "./util";
import { editNativeGasRegion, nativeGasRegions, type NativeGasRegion, type NativeGasValues } from "./native-gas-records.mjs";

/** Source-data controls shared by wind and string performances, without a browser gas model. */
export function createNativeAtmosphereControls(readSource: () => string, writeSource: (source:string) => void,
  draftChanged: () => void) {
  const selection = el("select", {"aria-label":"Gas region in the native performance"});
  const kind = el("select", {"aria-label":"Native gas model"},
    el("option", {value:"moist-air"}, "Moist air — relative humidity"),
    el("option", {value:"mixture"}, "Declared species mixture"),
    el("option", {value:"thermochemical"}, "Attached thermochemical species"));
  const field = (label:string) => el("input", {type:"text", inputmode:"decimal", spellcheck:"false", "aria-label":label});
  const temperature = field("Gas temperature in kelvin"), pressure = field("Absolute gas pressure in pascals");
  const humidity = field("Relative humidity as a fraction");
  const basis = el("select", {"aria-label":"Gas fraction basis"},
    el("option", {value:"mole"}, "Mole fractions"), el("option", {value:"mass"}, "Mass fractions"));
  const composition = el("div");
  const add = el("button", {type:"button", class:"btn"}, "Add species");
  const common = el("input", {type:"checkbox"});
  const commonLabel = el("label", {class:"small nw-decimation"}, common,
    " Apply this absolute pressure to every editable region");
  const apply = el("button", {type:"button", class:"btn primary"}, "Apply atmosphere");
  const cancel = el("button", {type:"button", class:"btn"}, "Discard changes");
  const status = el("p", {class:"small muted", role:"status", "aria-live":"polite"});
  const humidityLabel = el("label", {}, "Relative humidity (0–1)", humidity);
  const mixture = el("div", {}, el("label", {}, "Composition basis", basis), composition, add);
  const modelHelp = el("p", {class:"small muted"});
  const fields = el("fieldset", {class:"nw-gas-fields"},
    el("legend", {}, "Regional atmosphere"),
    el("div", {class:"nw-gas-policy"}, el("label", {}, "Region", selection), el("label", {}, "Gas model", kind)),
    el("div", {class:"nw-gas-policy"}, el("label", {}, "Temperature (K)", temperature),
      el("label", {}, "Absolute pressure (Pa)", pressure), humidityLabel),
    mixture, modelHelp, commonLabel,
    el("div", {class:"btn-row"}, apply, cancel), status);
  const empty = el("p", {class:"small muted"}, "No editable primitive gas records were found. Source-pack references, inherited aliases and spatial-field declarations remain in the performance editor.");
  const element = el("section", {class:"nw-gas-controls nw-atmosphere", "aria-label":"Edit physical atmosphere"},
    el("p", {class:"small muted"}, "Change a region's temperature, absolute pressure or composition, then render to hear the native physics. Geometry, excitation and pressure-to-PCM scale stay in the performance. Apply or discard these edits before continuing."), fields, empty);
  let regions: NativeGasRegion[] = [], selected: NativeGasRegion | undefined;
  let snapshot = "", locked = false, dirty = false;
  let rows: {species:HTMLInputElement; fraction:HTMLInputElement; element:HTMLElement}[] = [];
  const controls = (): void => {
    fields.disabled = locked;
    selection.disabled = locked || dirty;
    const legacy = kind.value === "moist-air";
    humidityLabel.hidden = !legacy;
    mixture.hidden = legacy;
    pressure.disabled = locked || (legacy && Boolean(selected?.inherited));
    common.disabled = locked || pressure.disabled;
    if (common.disabled) common.checked = false;
    add.disabled = locked || rows.length >= (kind.value === "mixture" ? 6 : 128);
    apply.disabled = locked || !dirty;
    cancel.disabled = locked || !dirty;
    modelHelp.textContent = legacy
      ? `The native moist-air owner derives the gas state from temperature, pressure and relative humidity.${selected?.inherited ? " This segment inherits ambient pressure; select ambient to change it." : ""}`
      : kind.value === "mixture"
        ? "Use the native species IDs N2, O2, Ar, CO2, H2O or He. Fractions must sum to 1; they are never normalized or converted between mole and mass automatically. Native material domains still apply."
        : "Use exact species IDs from the attached species packs. Keep their NASA-9 and transport associations below. Every selected species needs its original sources; the native model derives its state and rejects missing data.";
  };
  const changed = (): void => {
    dirty = true; status.textContent = "Atmosphere changes are pending. Apply them to the performance or discard them.";
    controls(); draftChanged();
  };
  const addRow = (value = {species:"", fraction:""}): void => {
    const species = el("input", {type:"text", spellcheck:"false", autocomplete:"off", "aria-label":"Exact gas species ID"});
    const fraction = field("Gas component fraction");
    species.value = value.species; fraction.value = value.fraction;
    const remove = el("button", {type:"button", class:"btn"}, "Remove species");
    const row = {species, fraction, element:el("div", {class:"nw-gas-species"},
      el("label", {}, "Species ID", species), el("label", {}, "Fraction (0–1)", fraction), remove)};
    species.addEventListener("input", changed); fraction.addEventListener("input", changed);
    remove.addEventListener("click", () => { rows = rows.filter(item => item !== row); row.element.remove(); changed(); });
    rows.push(row); composition.append(row.element);
  };
  const load = (): void => {
    selected = regions[Number(selection.value)];
    dirty = false; common.checked = false; rows = []; composition.replaceChildren();
    if (selected) {
      const value = selected.values;
      kind.value = value.kind; temperature.value = value.temperature; pressure.value = value.pressure;
      humidity.value = value.humidity ?? ""; basis.value = value.basis ?? "mole";
      value.components.forEach(addRow);
    }
    status.textContent = "Changes update the original gas records. Native boundary conditions may require one common static pressure.";
    controls();
  };
  const refresh = (source:string): void => {
    if (snapshot === source) return;
    const oldLine = selected?.line;
    snapshot = source; regions = nativeGasRegions(source);
    selection.replaceChildren(...regions.map((region, index) => el("option", {value:String(index)}, region.label)));
    selection.value = String(Math.max(0, regions.findIndex(region => region.line === oldLine)));
    fields.hidden = !regions.length; empty.hidden = Boolean(regions.length);
    load();
  };
  selection.addEventListener("change", load);
  kind.addEventListener("change", () => {
    if (kind.value !== "moist-air" && !rows.length) addRow();
    changed();
  });
  for (const input of [temperature, pressure, humidity]) input.addEventListener("input", changed);
  basis.addEventListener("change", changed); common.addEventListener("change", changed);
  add.addEventListener("click", () => { addRow(); changed(); });
  cancel.addEventListener("click", () => { load(); draftChanged(); });
  apply.addEventListener("click", () => {
    if (!selected || locked) return;
    const values: NativeGasValues = kind.value === "moist-air"
      ? {kind:"moist-air", temperature:temperature.value.trim(), pressure:pressure.value.trim(), humidity:humidity.value.trim(), components:[]}
      : {kind:kind.value as "mixture" | "thermochemical", temperature:temperature.value.trim(), pressure:pressure.value.trim(),
        basis:basis.value as "mole" | "mass", components:rows.map(row => ({species:row.species.value.trim(), fraction:row.fraction.value.trim()}))};
    try {
      const next = editNativeGasRegion(readSource(), selected, values, common.checked);
      dirty = false;
      writeSource(next);
      load();
      status.textContent = "Atmosphere applied to the native performance. Render again to hear it.";
      draftChanged();
    } catch (error) { status.textContent = error instanceof Error ? error.message : String(error); }
  });
  return {element, refresh, hasDraft:() => dirty, setDisabled(disabled:boolean) { locked = disabled; controls(); }};
}
