"use strict";
(() => {
  const $ = id => document.getElementById(id);
  const studio = window.COFdrawStudio;
  const dialog = $("solubility-dialog");
  const run = $("solubility-run"), status = $("solubility-status");
  const output = $("solubility-output"), downloads = $("solubility-exports");
  let smiles, revision = 0;

  function node(tag, content, className) {
    const el = document.createElement(tag);
    if (content !== undefined) el.textContent = content;
    if (className) el.className = className;
    return el;
  }
  function amount(value) {
    return value >= 0.01 && value < 10000
      ? Number(value.toPrecision(3)).toLocaleString(undefined, {maximumFractionDigits: 4})
      : value.toExponential(1);
  }
  function doiLink(doi, label) {
    const link = node("a", label || doi);
    link.href = `https://doi.org/${doi}`;
    link.target = "_blank";
    link.rel = "noopener noreferrer";
    return link;
  }
  function renderMeasurement(measurement) {
    const form = measurement.solid_form ? ` · ${measurement.solid_form}` : "";
    const method = measurement.method?.replaceAll("_", " ") || "equilibrium";
    const p = node("p", `${amount(measurement.mg_per_ml_solvent)} mg/mL at ${(measurement.temperature_k - 273.15).toFixed(1)} °C · ${method}${form} · `);
    p.append(doiLink(measurement.source_doi, "Study"));
    return p;
  }
  function renderSolvent(row) {
    const tr = node("tr");
    const name = node("th", row.solvent);
    name.scope = "row";
    tr.append(name);
    const origin = node("td");
    const disclosure = node("details");
    if (row.status === "reported") {
      disclosure.append(node("summary", "Reported", "solubility-reported"));
      for (const measurement of row.measurements) disclosure.append(renderMeasurement(measurement));
    } else {
      const prediction = row.evaluation;
      disclosure.append(node("summary", "Estimated", "solubility-estimated"));
      disclosure.append(node("p", `${row.analogs.length} measured analogs · ${prediction.predicted} held-out predictions · mean error ${prediction.mae_log10.toFixed(2)} log units`));
    }
    origin.append(disclosure);
    const concentration = node("td", amount(row.mg_per_ml_solvent));
    const mass = node("td", amount(row.mass_mg));
    if (row.screening_band) {
      const {lower, upper} = row.screening_band;
      concentration.append(node("small", `${amount(lower.mg_per_ml_solvent)}–${amount(upper.mg_per_ml_solvent)}`));
      mass.append(node("small", `${amount(lower.mass_mg)}–${amount(upper.mass_mg)}`));
    }
    tr.append(origin, concentration, mass);
    return tr;
  }
  function render(report) {
    const reported = report.solvents.filter(row => row.status === "reported");
    const estimated = report.solvents.filter(row => row.status === "estimated");
    const unavailable = report.solvents.filter(row => row.status === "unavailable");
    output.replaceChildren();
    if (reported.length || estimated.length) {
      const table = node("table", undefined, "solubility-table");
      const scroll = node("div", undefined, "solubility-table-scroll");
      const header = node("tr");
      for (const label of ["Solvent", "Evidence", "mg/mL solvent", `mg in ${amount(report.solvent_ml)} mL`]) {
        const th = node("th", label); th.scope = "col"; header.append(th);
      }
      const thead = node("thead"), body = node("tbody");
      thead.append(header);
      for (const row of [...reported, ...estimated]) body.append(renderSolvent(row));
      table.append(thead, body); scroll.append(table); output.append(scroll);
    }
    if (report.preparation_loadings?.length) {
      const section = node("section", undefined, "solubility-loadings");
      section.append(node("h4", "Reported preparation loads"));
      for (const load of report.preparation_loadings) {
        const p = node("p", `${load.solvent} · ${load.mass_mg} mg in ${load.solvent_ml} mL (${amount(load.mg_per_ml_solvent)} mg/mL) · `);
        p.append(doiLink(load.source_doi, load.source_locator));
        section.append(p);
      }
      section.append(node("p", "These are reported dissolution loads. Saturation and preparation temperature were unspecified.", "small-note"));
      output.append(section);
    }
    if (unavailable.length) {
      const details = node("details", undefined, "solubility-notes");
      details.open = !reported.length && !estimated.length;
      details.append(node("summary", `Unavailable solvents (${unavailable.length})`));
      const ul = node("ul");
      for (const row of unavailable) ul.append(node("li", `${row.solvent} · ${row.reason}`));
      details.append(ul); output.append(details);
    }
    const details = node("details", undefined, "solubility-notes");
    details.append(node("summary", "Method and assumptions"));
    const source = node("p", `${report.method.name} · `);
    source.append(doiLink(report.method.dataset_doi, report.method.dataset));
    details.append(source);
    for (const limitation of report.limitations) details.append(node("p", limitation));
    output.append(details);
    for (const [extension, path] of Object.entries(report.downloads || {})) {
      const button = node("button", `Save ${extension.toUpperCase()}`, "button outline");
      button.addEventListener("click", async () => {
        try { await studio.saveExport(path, extension); }
        catch (error) { status.textContent = error.message; }
      });
      downloads.append(button);
    }
    status.textContent = `${report.formula} · ${reported.length} reported · ${estimated.length} estimated`;
  }
  async function screen() {
    const token = ++revision;
    run.disabled = true; output.replaceChildren(); downloads.replaceChildren();
    status.textContent = "Screening solvents…";
    try {
      if (!smiles) throw new Error("Choose or draw a monomer first.");
      const result = await studio.api("science/solubility", {smiles, solvent_ml: $("solubility-volume").valueAsNumber});
      if (token === revision && dialog.open) render(result);
    } catch (error) {
      if (token === revision && dialog.open) status.textContent = error.message;
    } finally {
      if (token === revision) run.disabled = false;
    }
  }
  $("solubility-form").addEventListener("submit", event => { event.preventDefault(); screen(); });
  window.COFdrawSolubility = {
    open(monomer, name) {
      smiles = monomer;
      $("solubility-monomer").textContent = name;
      $("solubility-volume").value = "1";
      dialog.showModal();
      screen();
    },
  };
  $("descriptor-solubility").addEventListener("click", () => {
    const block = $("descriptor-block").value;
    window.COFdrawSolubility.open($(`smiles-${block}`).value.trim(), $(`label-${block}`).value || `Building block ${block.toUpperCase()}`);
  });
})();
