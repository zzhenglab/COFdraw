/* Small local 2D graph editor. Chemical validation belongs to the RDKit API. */
"use strict";
class MonomerSketcher {
  constructor(svg, onChange = () => {}, onError = () => {}, onSelect = () => {}) {
    this.svg = svg; this.onChange = onChange; this.onError = onError; this.graph = {atoms: [], bonds: []};
    this.onSelect = onSelect; this.tool = "bond"; this.element = "C"; this.order = 1; this.kind = null;
    this.selected = null; this.hover = null;
    this.history = []; this.future = []; this.drag = null; this.scale = 35; this.center = [0, 0];
    svg.addEventListener("pointerdown", event => { try { this.down(event); } catch (error) { this.onError(error.message); } });
    svg.addEventListener("pointermove", event => this.move(event));
    svg.addEventListener("pointerup", event => { try { this.up(event); } catch (error) { this.drag = null; this.draw(); this.onError(error.message); } });
    svg.addEventListener("pointercancel", () => this.cancelDrag());
    svg.addEventListener("pointerleave", () => { this.hover = null; if (!this.drag) this.draw(); });
    svg.addEventListener("keydown", event => {
      try { this.key(event); } catch (error) { this.onError(error.message); }
    });
    this.draw();
  }
  clone(value) { return JSON.parse(JSON.stringify(value)); }
  remember(previous = this.clone(this.graph)) { this.history.push(previous); if (this.history.length > 100) this.history.shift(); this.future = []; }
  changed() { this.draw(); this.onChange(this.exportGraph()); }
  exportGraph() { return this.clone(this.graph); }
  load(graph) {
    if (!graph || !Array.isArray(graph.atoms) || !Array.isArray(graph.bonds)) throw new Error("The editor needs atoms and bonds.");
    if (graph.atoms.length > 300) throw new Error("The monomer editor supports up to 300 atoms.");
    const atoms = graph.atoms.map(a => ({...a, element: a.element || a.symbol || "C", x: Number(a.x), y: Number(a.y), charge: Number(a.charge || 0)}));
    if (atoms.some(a => !Number.isFinite(a.x) || !Number.isFinite(a.y))) throw new Error("Atom coordinates must be finite.");
    const bonds = graph.bonds.map(b => ({...b, a: Number(b.a), b: Number(b.b), order: Number(b.order)}));
    if (bonds.some(b => !Number.isInteger(b.a) || !Number.isInteger(b.b) || b.a === b.b || !atoms[b.a] || !atoms[b.b] || ![1, 1.5, 2, 3].includes(b.order) || (b.kind && (b.kind !== "dative" || b.order !== 1)))) throw new Error("Invalid atom reference or bond order.");
    this.remember(); this.graph = {...this.clone(graph), atoms, bonds}; this.selected = null; this.hover = null; this.fit(); this.changed();
  }
  clear() { this.remember(); this.graph = {atoms: [], bonds: []}; this.selected = null; this.hover = null; this.fit(); this.changed(); }
  undo() { if (!this.history.length) return; this.future.push(this.clone(this.graph)); this.graph = this.history.pop(); this.selected = null; this.hover = null; this.fit(); this.changed(); }
  redo() { if (!this.future.length) return; this.history.push(this.clone(this.graph)); this.graph = this.future.pop(); this.selected = null; this.hover = null; this.fit(); this.changed(); }
  setTool(tool, value) {
    this.cancelDrag(); this.tool = tool;
    if (tool === "atom") this.element = value || "C";
    if (tool === "bond") { this.kind = value === "dative" ? value : null; this.order = this.kind ? 1 : Number(value || 1); }
    this.svg.dataset.tool = tool;
    this.svg.dispatchEvent(new CustomEvent("sketchtool", {detail: {tool, value: String(value || "")}}));
    this.draw();
  }
  key(event) {
    if (this.svg.getAttribute("aria-disabled") === "true") return;
    const key = event.key.toLowerCase();
    if (this.drag) { if (key === "escape") { event.preventDefault(); this.cancelDrag(); } return; }
    if (event.ctrlKey || event.metaKey) {
      if (key === "z" || key === "y") {
        event.preventDefault();
        if (key === "y" || event.shiftKey) this.redo();
        else this.undo();
      }
      return;
    }
    if (event.altKey) return;
    if (key === "delete" || key === "backspace") { event.preventDefault(); this.removeSelected(); }
    else if (key === "escape") { this.selected = null; this.draw(); }
    else if (["1", "2", "3", "d"].includes(key)) {
      event.preventDefault();
      const target = this.hover || this.selected;
      this.setTool("bond", key === "d" ? "dative" : key);
      if (target?.kind === "bond") this.changeBond(target.index);
    } else if (["c", "n", "o", "s", "f", "p", "b", "h", "i"].includes(key)) {
      event.preventDefault();
      const target = this.hover || this.selected;
      this.setTool("atom", key.toUpperCase());
      if (target?.kind === "atom") this.placeElement(this.graph.atoms[target.index], target.index);
    }
  }
  cancelDrag() {
    if (!this.drag) return;
    if (this.drag.previous) this.graph = this.drag.previous;
    this.drag = null; this.draw();
  }
  fit() {
    if (!this.graph.atoms.length) { this.scale = 35; this.center = [0, 0]; return; }
    const xs = this.graph.atoms.map(a => a.x), ys = this.graph.atoms.map(a => a.y);
    const width = Math.max(...xs) - Math.min(...xs), height = Math.max(...ys) - Math.min(...ys);
    this.center = [(Math.min(...xs) + Math.max(...xs)) / 2, (Math.min(...ys) + Math.max(...ys)) / 2];
    this.scale = Math.min(38, 540 / Math.max(3, width), 270 / Math.max(3, height));
  }
  screen(atom) { return [320 + (atom.x - this.center[0]) * this.scale, 180 - (atom.y - this.center[1]) * this.scale]; }
  position(event) {
    const p = new DOMPoint(event.clientX, event.clientY).matrixTransform(this.svg.getScreenCTM().inverse());
    return {x: (p.x - 320) / this.scale + this.center[0], y: (180 - p.y) / this.scale + this.center[1]};
  }
  atomAt(p) { let found = -1, distance = 13 / this.scale; this.graph.atoms.forEach((a, i) => { const d = Math.hypot(a.x - p.x, a.y - p.y); if (d < distance) { found = i; distance = d; } }); return found; }
  bondAt(p) {
    let found = -1, distance = 8 / this.scale;
    this.graph.bonds.forEach((b, i) => { const a = this.graph.atoms[b.a], c = this.graph.atoms[b.b]; const dx = c.x - a.x, dy = c.y - a.y, den = dx * dx + dy * dy; if (!den) return; const t = Math.max(0, Math.min(1, ((p.x - a.x) * dx + (p.y - a.y) * dy) / den)); const d = Math.hypot(p.x - a.x - t * dx, p.y - a.y - t * dy); if (d < distance) { found = i; distance = d; } }); return found;
  }
  hasStereo() { return this.graph.atoms.some(a => a.chiral_tag && !["CHI_UNSPECIFIED", "UNSPECIFIED", "0"].includes(String(a.chiral_tag))) || this.graph.bonds.some(b => (b.stereo && !["STEREONONE", "NONE", "0"].includes(String(b.stereo))) || (b.direction && !["NONE", "0"].includes(String(b.direction)))); }
  connectivityEditable() { if (this.hasStereo()) throw new Error("This monomer has assigned stereochemistry. Move atoms here, or edit its stereospecific SMILES / exported MOL in a full chemical editor to change connectivity."); }
  addAtom(p, element = "C") { if (this.graph.atoms.length >= 300) throw new Error("The editor limit is 300 atoms."); this.graph.atoms.push({element, x: p.x, y: p.y, charge: 0}); return this.graph.atoms.length - 1; }
  addBond(a, b, order = this.order, kind = null) {
    if (a === b) return;
    const old = this.graph.bonds.findIndex(x => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    const bond = {a, b, order, ...(kind ? {kind} : {})};
    if (old >= 0) this.graph.bonds[old] = bond;
    else this.graph.bonds.push(bond);
  }
  changeBond(index) {
    this.connectivityEditable();
    const bond = this.graph.bonds[index];
    if (!bond || (!this.kind && !bond.kind && bond.order === this.order)) return;
    this.remember();
    // Reapplying the arrow tool reverses an existing coordination bond.
    const [a, b] = this.kind && bond.kind ? [bond.b, bond.a] : [bond.a, bond.b];
    this.addBond(a, b, this.order, this.kind);
    this.selected = {kind: "bond", index};
    this.changed();
  }
  placeElement(p, index) {
    this.connectivityEditable();
    if (index < 0 && this.graph.atoms.length >= 300) throw new Error("The editor limit is 300 atoms.");
    if (this.graph.atoms[index]?.element === this.element) { this.selected = {kind: "atom", index}; this.draw(); return; }
    this.remember();
    if (index < 0) index = this.addAtom(p, this.element);
    else {
      const {x, y, charge} = this.graph.atoms[index];
      this.graph.atoms[index] = {element: this.element, x, y, charge};
    }
    this.selected = {kind: "atom", index};
    this.changed();
  }
  deleteAtom(index) { this.graph.atoms.splice(index, 1); this.graph.bonds = this.graph.bonds.filter(b => b.a !== index && b.b !== index).map(b => ({...b, a: b.a > index ? b.a - 1 : b.a, b: b.b > index ? b.b - 1 : b.b})); }
  removeSelected() { if (!this.selected) return; this.connectivityEditable(); this.remember(); if (this.selected.kind === "atom") this.deleteAtom(this.selected.index); else this.graph.bonds.splice(this.selected.index, 1); this.selected = null; this.hover = null; this.changed(); }
  charge(value) { if (this.selected?.kind !== "atom") return false; this.remember(); this.graph.atoms[this.selected.index].charge = Math.max(-4, Math.min(8, Number(value))); this.changed(); return true; }
  ring(p, existing, count, aromatic) {
    if (this.graph.atoms.length + count > 300) return;
    this.remember(); const radius = 1.5 / (2 * Math.sin(Math.PI / count));
    const anchor = existing >= 0 ? this.graph.atoms[existing] : p;
    const angle = existing >= 0 ? this.openAngle(existing) : 0;
    const center = existing >= 0 ? this.endAt(anchor, angle, radius) : p;
    const ids = [];
    for (let i = 0; i < count; i++) ids.push(i === 0 && existing >= 0 ? existing : this.addAtom(this.endAt(center, angle + Math.PI - i * 2 * Math.PI / count, radius)));
    for (let i = 0; i < count; i++) this.addBond(ids[i], ids[(i + 1) % count], aromatic && i % 2 === 0 ? 2 : 1);
    this.selected = {kind: "atom", index: ids[0]}; this.changed();
  }
  openAngle(index) {
    if (index < 0) return Math.PI / 6;
    const atom = this.graph.atoms[index];
    const neighbors = this.graph.bonds.filter(b => b.a === index || b.b === index).map(b => {
      const p = this.graph.atoms[b.a === index ? b.b : b.a];
      return {angle: Math.atan2(p.y - atom.y, p.x - atom.x), order: b.order};
    });
    if (!neighbors.length) return Math.PI / 6;
    const [first] = neighbors;
    if (neighbors.length === 1 && (Math.max(this.order, first.order) === 3 || (this.order === 2 && first.order === 2))) return first.angle + Math.PI;
    const directions = neighbors.length === 1
      ? [first.angle + 2 * Math.PI / 3, first.angle - 2 * Math.PI / 3]
      : Array.from({length: 12}, (_, i) => first.angle + i * Math.PI / 6);
    const score = angle => {
      const p = this.endAt(atom, angle);
      const room = Math.min(1.5, ...this.graph.atoms.map(a => Math.hypot(a.x - p.x, a.y - p.y)));
      const gap = Math.min(...neighbors.map(a => Math.acos(Math.max(-1, Math.min(1, Math.cos(a.angle - angle))))));
      return gap + room;
    };
    return directions.reduce((best, candidate) => score(candidate) > score(best) + 1e-8 ? candidate : best);
  }
  endAt(start, angle, length = 1.5) {
    return {x: start.x + length * Math.cos(angle), y: start.y + length * Math.sin(angle)};
  }
  hit(p) {
    const atom = this.atomAt(p);
    if (atom >= 0) return {kind: "atom", index: atom};
    const bond = this.bondAt(p);
    return bond >= 0 ? {kind: "bond", index: bond} : null;
  }
  down(event) {
    if (event.button !== 0) return; event.preventDefault(); this.svg.focus();
    const p = this.position(event), atom = this.atomAt(p), bond = atom < 0 ? this.bondAt(p) : -1;
    if (this.tool !== "select") this.connectivityEditable();
    if (this.tool.startsWith("ring")) { this.ring(p, atom, this.tool === "ring5" ? 5 : 6, this.tool === "ring-aromatic"); return; }
    if (this.tool === "erase") { this.selected = this.hit(p); this.removeSelected(); return; }
    if (this.tool === "atom") { this.placeElement(p, atom); return; }
    if (this.tool === "select") {
      this.selected = this.hit(p);
      if (atom >= 0) {
        const a = this.graph.atoms[atom];
        this.drag = {kind: "move", atom, dx: a.x - p.x, dy: a.y - p.y, previous: this.exportGraph()};
        this.svg.setPointerCapture(event.pointerId);
      }
      this.draw();
      return;
    }
    if (bond >= 0) { this.changeBond(bond); return; }
    const start = atom >= 0 ? {...this.graph.atoms[atom]} : p;
    this.drag = {kind: "bond", atom, start, origin: p, end: this.endAt(start, this.openAngle(atom)), target: -1};
    this.svg.setPointerCapture(event.pointerId); this.draw();
  }
  move(event) {
    if (this.svg.getAttribute("aria-disabled") === "true") return;
    const p = this.position(event), drag = this.drag;
    if (!drag) {
      const hit = this.hit(p);
      if (hit?.kind !== this.hover?.kind || hit?.index !== this.hover?.index) { this.hover = hit; this.draw(); }
      return;
    }
    if (drag.kind === "move") {
      const atom = this.graph.atoms[drag.atom];
      atom.x = p.x + drag.dx; atom.y = p.y + drag.dy; this.draw(); return;
    }
    const target = this.atomAt(p), length = Math.hypot(p.x - drag.start.x, p.y - drag.start.y);
    drag.target = target === drag.atom ? -1 : target;
    if (drag.target >= 0) drag.end = {...this.graph.atoms[drag.target]};
    else if (Math.hypot(p.x - drag.origin.x, p.y - drag.origin.y) * this.scale < 5) drag.end = this.endAt(drag.start, this.openAngle(drag.atom));
    else {
      const angle = Math.atan2(p.y - drag.start.y, p.x - drag.start.x);
      const step = Math.PI / 6;
      drag.end = event.altKey ? this.endAt(drag.start, angle, Math.max(.5, length)) : this.endAt(drag.start, Math.round(angle / step) * step);
    }
    this.draw();
  }
  up(event) {
    if (!this.drag) return; const drag = this.drag;
    if (this.svg.hasPointerCapture(event.pointerId)) this.svg.releasePointerCapture(event.pointerId);
    this.drag = null;
    if (drag.kind === "move") {
      const before = drag.previous.atoms[drag.atom], after = this.graph.atoms[drag.atom];
      if (Math.hypot(before.x - after.x, before.y - after.y) < .001) { this.graph = drag.previous; this.draw(); return; }
      this.remember(drag.previous); this.changed(); return;
    }
    if (this.graph.atoms.length + Number(drag.atom < 0) + Number(drag.target < 0) > 300) throw new Error("The editor limit is 300 atoms.");
    this.remember();
    const a = drag.atom >= 0 ? drag.atom : this.addAtom(drag.start);
    const b = drag.target >= 0 ? drag.target : this.addAtom(drag.end);
    this.addBond(a, b, this.order, this.kind); this.selected = {kind: "atom", index: b}; this.changed();
  }
  atomLabel(atom) {
    const hydrogens = atom.explicit_hs || 0;
    const charge = atom.charge ? (Math.abs(atom.charge) > 1 ? Math.abs(atom.charge) : "") + (atom.charge > 0 ? "+" : "−") : "";
    return `${atom.isotope || ""}${atom.element}${hydrogens ? "H" + (hydrogens > 1 ? hydrogens : "") : ""}${charge}`;
  }
  node(name, attrs = {}, text) { const e = document.createElementNS("http://www.w3.org/2000/svg", name); Object.entries(attrs).forEach(([k, v]) => e.setAttribute(k, v)); if (text !== undefined) e.textContent = text; return e; }
  draw() {
    const s = this.svg; s.replaceChildren();
    if (!this.graph.atoms.length && !this.drag) { s.append(this.node("text", {x:320,y:165,"text-anchor":"middle",class:"sketch-placeholder"}, "Choose a ring or draw a bond")); s.append(this.node("text", {x:320,y:188,"text-anchor":"middle",class:"sketch-placeholder small"}, "Click an atom to change its element. Drag to connect.")); }
    this.graph.bonds.forEach((b, index) => {
      const a = this.screen(this.graph.atoms[b.a]), c = this.screen(this.graph.atoms[b.b]); const d = Math.hypot(c[0]-a[0],c[1]-a[1]) || 1, nx = -(c[1]-a[1])/d, ny = (c[0]-a[0])/d;
      if (this.selected?.kind === "bond" && this.selected.index === index) s.append(this.node("line",{x1:a[0],y1:a[1],x2:c[0],y2:c[1],stroke:"#b8e2d5","stroke-width":12,"stroke-linecap":"round"}));
      else if (this.hover?.kind === "bond" && this.hover.index === index) s.append(this.node("line",{x1:a[0],y1:a[1],x2:c[0],y2:c[1],class:"sketch-bond-hover"}));
      if (b.kind === "dative") {
        const ux = ny, uy = -nx, acceptor = this.graph.atoms[b.b];
        const label = this.atomLabel(acceptor);
        const inset = label === "C" ? 1 : Math.min((4 * label.length + 2) / Math.abs(ux), 12 / Math.abs(uy));
        const x = c[0] - ux * Math.min(inset, d * .45), y = c[1] - uy * Math.min(inset, d * .45);
        s.append(this.node("line",{x1:a[0],y1:a[1],x2:x,y2:y,class:"sketch-dative"}));
        s.append(this.node("path",{d:`M ${x-ux*7+nx*3} ${y-uy*7+ny*3} L ${x} ${y} L ${x-ux*7-nx*3} ${y-uy*7-ny*3}`,class:"sketch-dative-arrow"}));
        return;
      }
      const offsets = b.order === 3 ? [-4,0,4] : b.order === 2 || b.order === 1.5 ? [-2.5,2.5] : [0];
      offsets.forEach((o,j) => s.append(this.node("line",{x1:a[0]+nx*o,y1:a[1]+ny*o,x2:c[0]+nx*o,y2:c[1]+ny*o,stroke:"#243c3a","stroke-width":1.8,"stroke-linecap":"round",...(b.order===1.5&&j===1?{"stroke-dasharray":"4 3"}:{})})));
    });
    this.graph.atoms.forEach((atom,index) => {
      const [x,y] = this.screen(atom), selected = this.selected?.kind === "atom" && this.selected.index === index;
      if (selected) s.append(this.node("circle",{cx:x,cy:y,r:12,fill:"#d9f1e8",stroke:"#2b8770"}));
      else if (this.hover?.kind === "atom" && this.hover.index === index) s.append(this.node("circle",{cx:x,cy:y,r:12,class:"sketch-atom-hover"}));
      const degree = this.graph.bonds.filter(b=>b.a===index||b.b===index).length;
      if (atom.element !== "C" || atom.charge || atom.isotope || atom.explicit_hs || !degree) {
        const label = this.atomLabel(atom);
        s.append(this.node("rect",{x:x-8*label.length/2,y:y-10,width:8*label.length,height:20,fill:"white"}));
        s.append(this.node("text",{x,y:y+5,"text-anchor":"middle",fill:({N:"#3564a6",O:"#b84444",S:"#9a760d",F:"#39834a",Cl:"#39834a",Br:"#934533"})[atom.element]||"#243c3a",class:"sketch-atom"},label));
      } else if (selected) s.append(this.node("circle",{cx:x,cy:y,r:2.5,fill:"#206b56"}));
    });
    if (this.drag?.kind === "bond") { const a=this.screen(this.drag.start),b=this.screen(this.drag.end); s.append(this.node("line",{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:"#31a483","stroke-width":2,"stroke-dasharray":"5 4"})); }
    this.onSelect(this.selected?.kind === "atom" ? this.graph.atoms[this.selected.index] : null);
  }
}
window.MonomerSketcher = MonomerSketcher;

function mountElementPicker(dialog, choose) {
  // Columns follow conventional periodic groups for the common elements shown.
  const elements = [
    [1, 1, 1, "H", "Hydrogen"],
    [2, 1, 3, "Li", "Lithium"], [2, 2, 4, "Be", "Beryllium"],
    [2, 13, 5, "B", "Boron"], [2, 14, 6, "C", "Carbon"], [2, 15, 7, "N", "Nitrogen"],
    [2, 16, 8, "O", "Oxygen"], [2, 17, 9, "F", "Fluorine"],
    [3, 1, 11, "Na", "Sodium"], [3, 2, 12, "Mg", "Magnesium"],
    [3, 13, 13, "Al", "Aluminium"], [3, 14, 14, "Si", "Silicon"], [3, 15, 15, "P", "Phosphorus"],
    [3, 16, 16, "S", "Sulfur"], [3, 17, 17, "Cl", "Chlorine"],
    [4, 1, 19, "K", "Potassium"], [4, 2, 20, "Ca", "Calcium"], [4, 3, 21, "Sc", "Scandium"],
    [4, 4, 22, "Ti", "Titanium"], [4, 5, 23, "V", "Vanadium"], [4, 6, 24, "Cr", "Chromium"],
    [4, 7, 25, "Mn", "Manganese"], [4, 8, 26, "Fe", "Iron"], [4, 9, 27, "Co", "Cobalt"],
    [4, 10, 28, "Ni", "Nickel"], [4, 11, 29, "Cu", "Copper"], [4, 12, 30, "Zn", "Zinc"],
    [4, 13, 31, "Ga", "Gallium"], [4, 14, 32, "Ge", "Germanium"], [4, 15, 33, "As", "Arsenic"],
    [4, 16, 34, "Se", "Selenium"], [4, 17, 35, "Br", "Bromine"],
    [5, 1, 37, "Rb", "Rubidium"], [5, 2, 38, "Sr", "Strontium"], [5, 3, 39, "Y", "Yttrium"],
    [5, 4, 40, "Zr", "Zirconium"], [5, 5, 41, "Nb", "Niobium"], [5, 6, 42, "Mo", "Molybdenum"],
    [5, 8, 44, "Ru", "Ruthenium"], [5, 9, 45, "Rh", "Rhodium"], [5, 10, 46, "Pd", "Palladium"],
    [5, 11, 47, "Ag", "Silver"], [5, 12, 48, "Cd", "Cadmium"], [5, 13, 49, "In", "Indium"],
    [5, 14, 50, "Sn", "Tin"], [5, 15, 51, "Sb", "Antimony"], [5, 16, 52, "Te", "Tellurium"], [5, 17, 53, "I", "Iodine"],
    [6, 1, 55, "Cs", "Caesium"], [6, 2, 56, "Ba", "Barium"], [6, 3, 57, "La", "Lanthanum"],
    [6, 4, 72, "Hf", "Hafnium"], [6, 5, 73, "Ta", "Tantalum"], [6, 6, 74, "W", "Tungsten"],
    [6, 7, 75, "Re", "Rhenium"], [6, 8, 76, "Os", "Osmium"], [6, 9, 77, "Ir", "Iridium"],
    [6, 10, 78, "Pt", "Platinum"], [6, 11, 79, "Au", "Gold"], [6, 12, 80, "Hg", "Mercury"],
    [6, 14, 82, "Pb", "Lead"], [6, 15, 83, "Bi", "Bismuth"],
  ];
  const grid = dialog.querySelector("#element-grid"), search = dialog.querySelector("#element-search");
  for (const [period, group, number, symbol, name] of elements) {
    const button = document.createElement("button");
    button.type = "button";
    button.style.gridRow = period; button.style.gridColumn = group;
    button.className = group >= 3 && group <= 12 ? "element-cell transition" : "element-cell";
    button.title = `${name} · atomic number ${number}`;
    button.setAttribute("aria-label", `${name} (${symbol})`);
    button.dataset.search = `${symbol} ${name} ${number}`.toLowerCase();
    const n = document.createElement("small"); n.textContent = number;
    const s = document.createElement("strong"); s.textContent = symbol;
    button.append(n, s);
    button.addEventListener("click", () => { dialog.close(); choose(symbol); });
    grid.append(button);
  }
  function filter() {
    const query = search.value.trim().toLowerCase().replace("aluminum", "aluminium").replace("cesium", "caesium");
    const matches = [...grid.children].filter(button => {
      const found = button.dataset.search.split(" ").some(term => term.startsWith(query));
      button.disabled = !found;
      return found;
    });
    dialog.querySelector("#element-match").textContent = matches.length ? "" : "No matching element in this picker.";
    return matches;
  }
  search.addEventListener("input", filter);
  search.addEventListener("keydown", event => {
    if (event.key !== "Enter") return;
    event.preventDefault();
    const matches = filter();
    if (matches.length === 1) matches[0].click();
    else if (matches.length) matches[0].focus();
  });
  return () => { search.value = ""; filter(); dialog.showModal(); search.focus(); };
}
window.mountElementPicker = mountElementPicker;
