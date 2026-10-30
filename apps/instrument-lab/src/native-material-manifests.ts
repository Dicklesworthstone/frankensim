import { el } from "./util";
import type { NativeMaterialManifest } from "./native-wind-client.mjs";

/** Source associations only: the Rust compiler owns every material property and identity. */
export function createNativeMaterialManifests(changed: () => void) {
  const list = el("div", {class: "nw-gas-components"});
  const add = el("button", {type: "button", class: "btn"}, "Add material manifest");
  const availability = el("p", {class: "small muted"});
  const element = el("fieldset", {class: "nw-gas-fields"},
    el("legend", {}, "Compile original material sources"),
    el("p", {class: "small muted"}, "Attach each original manifest and the source files it names, preserving their relative paths. Choose the pack path used by the performance or gas association. Remove any precompiled attachment at that output path. Edit the attached source tables to change physical data; rendering recompiles them through the native material compiler."),
    availability, list, add);
  const identityKeys = ["chemistry", "phase", "process", "revision"] as const;
  type Row = {element: HTMLElement; pack: HTMLInputElement; manifest: HTMLSelectElement;
    named: HTMLInputElement; state: Record<typeof identityKeys[number], HTMLInputElement>};
  let rows: Row[] = [];
  let paths: readonly string[] = [];
  const choices = (input: HTMLSelectElement, selected = input.value): void => {
    input.replaceChildren(el("option", {value: ""}, "Choose an attached manifest…"),
      ...paths.map(path => el("option", {value: path}, path)));
    if (selected && !paths.includes(selected)) input.append(el("option", {value: selected}, `Missing: ${selected}`));
    input.value = selected;
  };
  const addRow = (value?: NativeMaterialManifest): void => {
    const pack = el("input", {type: "text", spellcheck: "false", "aria-label": "Compiled material pack path"});
    pack.value = value?.packPath ?? "";
    const manifest = el("select", {"aria-label": "Original attached material manifest"});
    choices(manifest, value?.manifestPath ?? "");
    const named = el("input", {type: "checkbox"});
    named.checked = value?.materialState !== undefined;
    const state = Object.fromEntries(identityKeys.map(key => {
      const input = el("input", {type: "text", spellcheck: "false", "aria-label": `Material ${key}`});
      input.value = value?.materialState ? String(value.materialState[key]) : "";
      return [key, input];
    })) as Record<typeof identityKeys[number], HTMLInputElement>;
    state.revision.inputMode = "numeric";
    const stateFields = el("fieldset", {class: "nw-gas-component"},
      el("legend", {}, "Explicit material identity"),
      ...identityKeys.map(key => el("label", {}, key, state[key])));
    const refreshState = (): void => { stateFields.hidden = stateFields.disabled = !named.checked; };
    const remove = el("button", {type: "button", class: "btn"}, "Remove binding");
    const row: Row = {pack, manifest, named, state, element: el("div", {class: "nw-gas-component"},
      el("label", {}, "Original manifest", manifest), el("label", {}, "Output pack path", pack),
      el("label", {}, named, " Bind claims to a named material card"),
      el("p", {class: "small muted"}, "Use the explicit chemistry, phase, process and revision required by your material source. Leave this off for original claims, model, interface or species packs. No identity or coefficient is inferred."),
      stateFields, remove)};
    named.addEventListener("change", () => { refreshState(); changed(); });
    manifest.addEventListener("change", changed);
    [pack, ...Object.values(state)].forEach(input => input.addEventListener("input", changed));
    remove.addEventListener("click", () => {
      rows = rows.filter(item => item !== row); row.element.remove();
      add.disabled = rows.length >= 32; changed();
    });
    refreshState(); rows.push(row); list.append(row.element); add.disabled = rows.length >= 32;
  };
  add.addEventListener("click", () => { if (rows.length < 32) { addRow(); changed(); } });
  return {
    element,
    read(): NativeMaterialManifest[] | undefined {
      if (!rows.length) return undefined;
      return rows.map(row => {
        if (!paths.includes(row.manifest.value) || !row.pack.value)
          throw new Error("Choose an attached manifest and the output pack path for each material binding.");
        const binding: NativeMaterialManifest = {packPath:row.pack.value, manifestPath:row.manifest.value};
        if (row.named.checked) {
          const {chemistry, phase, process, revision} = row.state;
          const number = Number(revision.value);
          if ([chemistry, phase, process].some(input => !input.value.trim())
              || !/^[0-9]+$/.test(revision.value) || !Number.isInteger(number) || number > 4294967295)
            throw new Error("Enter the explicit chemistry, phase, process and unsigned revision for each named material card.");
          binding.materialState = {chemistry:chemistry.value, phase:phase.value, process:process.value, revision:number};
        }
        return binding;
      });
    },
    packPaths: (): string[] => rows.map(row => row.pack.value).filter(Boolean),
    identity: (): string => JSON.stringify(rows.map(row => [row.pack.value, row.manifest.value,
      row.named.checked, ...identityKeys.map(key => row.state[key].value)])),
    setAttachments(next: readonly string[]): void { paths = next; rows.forEach(row => choices(row.manifest)); },
    setAvailable(available: boolean): void {
      availability.textContent = available
        ? "The local material compiler is available. Original source inputs and bindings are retained in saved experiments."
        : "To render these bindings, configure the local server with --material-compiler /absolute/path/xtask. Source editing and experiment saving remain available.";
    },
    setDisabled(disabled: boolean): void { element.disabled = disabled; },
    reset(): void { rows = []; list.replaceChildren(); add.disabled = false; },
    load(bindings?: readonly NativeMaterialManifest[]): void {
      rows = []; list.replaceChildren(); add.disabled = false;
      bindings?.forEach(binding => addRow(binding));
    },
  };
}
