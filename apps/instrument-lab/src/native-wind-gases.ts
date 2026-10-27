import { el } from "./util";
import type { NativeWindGasComponent, NativeWindGasSources } from "./native-wind-client.mjs";

interface GasControls {
  element: HTMLElement;
  read(): NativeWindGasSources | undefined;
  identity(): string;
  setAttachments(paths: readonly string[]): void;
  setDisabled(disabled: boolean): void;
  reset(): void;
}

/** Explicit source associations only; the native material model resolves every coefficient. */
export function createNativeGasControls(changed: () => void): GasControls {
  const enabled = el("input", { type: "checkbox", id: "native-wind-gas-sources" });
  const list = el("div", { class: "nw-gas-components" });
  const add = el("button", { type: "button", class: "btn" }, "Add gas component");
  const phase = el("select", { "aria-label": "Gas phase policy" },
    el("option", { value: "" }, "Choose phase policy…"),
    el("option", { value: "frozen" }, "Frozen gas composition"),
    el("option", { value: "liquid-water" }, "Check liquid-water saturation"));
  const water = el("input", { type: "text", spellcheck: "false", autocomplete: "off",
    placeholder: "Species ID, or none", "aria-label": "Water species ID or none" });
  const fields = el("fieldset", { class: "nw-gas-fields" },
    el("legend", {}, "Thermochemical gas sources"),
    el("p", { class: "small muted" }, "Attach each species pack, its NASA-9 model pack, and its pure transport claims pack. Select their exact paths below. The performance keeps the species IDs, mixture fractions, temperature and pressure in its thermochemical gas records."),
    list, add,
    el("div", { class: "nw-gas-policy" },
      el("label", {}, "Phase policy", phase),
      el("label", {}, "Water identity", water)),
    el("p", { class: "small muted" }, "Use the original native transport claims pack, not a wrapped material-card pack. Missing coefficients, incompatible source IDs and states outside the source domains are refused by the native model."));
  const element = el("section", { class: "nw-gas-controls", "aria-label": "Gas material sources" },
    el("label", { for: "native-wind-gas-sources", class: "small nw-decimation" }, enabled,
      " Use attached thermochemical gas sources"), fields);
  const keys = ["speciesPath", "nasa9Path", "transportPath"] as const;
  const labels = ["Species pack", "NASA-9 model pack", "Pure transport claims pack"] as const;
  type Row = { element: HTMLElement; inputs: Record<typeof keys[number], HTMLSelectElement> };
  let rows: Row[] = [];
  let paths: readonly string[] = [];
  let locked = false;

  const choices = (input: HTMLSelectElement): void => {
    const selected = input.value;
    input.replaceChildren(el("option", { value: "" }, "Select an attached file…"),
      ...paths.map(name => el("option", { value: name }, name)));
    // A renamed/removed source must remain visibly unresolved, never silently rebound.
    if (selected && !paths.includes(selected)) input.append(el("option", { value: selected }, `Missing: ${selected}`));
    input.value = selected;
  };
  const controls = (): void => {
    enabled.disabled = locked;
    fields.hidden = !enabled.checked;
    fields.disabled = locked || !enabled.checked;
    add.disabled = locked || rows.length >= 128;
  };
  const components = (): NativeWindGasComponent[] => rows.map(row => ({
    speciesPath: row.inputs.speciesPath.value,
    nasa9Path: row.inputs.nasa9Path.value,
    transportPath: row.inputs.transportPath.value,
  }));
  const addRow = (): void => {
    if (rows.length >= 128) return;
    const inputs = { speciesPath: el("select"), nasa9Path: el("select"), transportPath: el("select") };
    const remove = el("button", { type: "button", class: "btn" }, "Remove component");
    const row: Row = { inputs, element: el("div", { class: "nw-gas-component" }) };
    keys.forEach((key, index) => {
      choices(inputs[key]);
      inputs[key].addEventListener("change", changed);
      row.element.append(el("label", {}, labels[index]!, inputs[key]));
    });
    remove.addEventListener("click", () => {
      rows = rows.filter(candidate => candidate !== row);
      row.element.remove();
      controls(); changed();
    });
    row.element.append(remove);
    rows.push(row); list.append(row.element); controls();
  };
  add.addEventListener("click", () => { addRow(); changed(); });
  enabled.addEventListener("change", () => { controls(); changed(); });
  phase.addEventListener("change", changed);
  water.addEventListener("input", changed);
  addRow();

  return {
    element,
    read() {
      if (!enabled.checked) return undefined;
      const sources = components();
      if (!sources.length || sources.some(component => keys.some(key => !component[key] || !paths.includes(component[key])))) {
        throw new Error("Choose all three attached source files for every gas component.");
      }
      if (phase.value !== "frozen" && phase.value !== "liquid-water") {
        throw new Error("Choose an explicit gas phase policy.");
      }
      const waterId = water.value.trim();
      if (!waterId) throw new Error("Enter the canonical water species ID, or none when water is absent.");
      const waterSpecies = waterId === "none" ? null : waterId;
      if (phase.value === "liquid-water" && waterSpecies === null) {
        throw new Error("The liquid-water policy needs the water species ID from your source pack.");
      }
      return {components:sources, phase:phase.value, waterSpecies};
    },
    identity: () => enabled.checked
      ? JSON.stringify({components:components(), phase:phase.value, water:water.value}) : "",
    setAttachments(next) {
      paths = [...new Set(next)];
      for (const row of rows) for (const key of keys) choices(row.inputs[key]);
    },
    setDisabled(disabled) { locked = disabled; controls(); },
    reset() {
      enabled.checked = false; phase.value = ""; water.value = "";
      rows = []; list.replaceChildren(); addRow(); controls(); changed();
    },
  };
}
