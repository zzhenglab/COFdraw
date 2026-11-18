# COFdraw agent entrypoint

The portable, installable skill is **[skills/cofdraw/SKILL.md](skills/cofdraw/SKILL.md)**. Read it when a user asks to use this GitHub repository to query local COF information, preserve source layouts, export chemical drawings, inspect periodic CIF/contact networks, calculate ideal PXRD or geometric bounds, screen finite-monomer SA Scores, run qualified local adsorption jobs, fit supplied isotherms, or retrieve saved reports. Its helper supports both a source checkout and the complete Windows release.

Install the `skills/cofdraw` folder, rather than this repository root, in hosts that support skills. A repository link does not itself install the application or grant an agent local execution tools.

Example request:

> Use https://github.com/zzhenglab/COFdraw and its skills/cofdraw/SKILL.md to draw the catalog TpPa pore in Standard style. Save PNG, SVG, editable chemical exports, and the provenance JSON.
