"""Small, reproducible comparisons of user-supplied periodic COF models."""
from __future__ import annotations

import csv
import io
import json
import re


def screen_structures(entries, *, include_porosity=False, probe_radius=1.82, samples=4000):
    from .cofdraw_characterization import structure_summary
    from .cofdraw_network import inspect_periodic_network
    from .cofdraw_science import pore_analysis

    if not isinstance(entries, list) or not 1 <= len(entries) <= 8:
        raise ValueError('Compare between one and eight structures')
    if len(json.dumps(entries, allow_nan=False)) > 6_000_000:
        raise ValueError('Batch input limit is 6 MB')
    results = []
    table = io.StringIO(newline='')
    fields = ['name', 'status', 'formula', 'density_g_cm3', 'cell_volume_angstrom3',
              'atom_count', 'inferred_network_dimensions', 'void_fraction',
              'geometric_void_volume_cm3_g', 'local_probe_exposed_surface_m2_g',
              'largest_sampled_included_sphere_diameter_angstrom', 'error']
    writer = csv.DictWriter(table, fieldnames=fields)
    writer.writeheader()
    for index, entry in enumerate(entries):
        name = str(entry.get('name', f'Structure {index+1}'))[:100]
        # Spreadsheet exports must not execute a user-provided formula.
        csv_name = "'" + name if re.match(r'^[\s\ufeff]*[=+@-]',name) or name.startswith(('\t','\r','\n')) else name
        row = {'name': name, 'status': 'error'}
        flat = {'name': csv_name, 'status': 'error'}
        try:
            summary = structure_summary(entry['structure'])
            network = inspect_periodic_network(entry['structure'])
            dimensions = sorted({c['periodic_rank'] for c in network['components']})
            row.update(status='ok', summary=summary, network=network,
                       inferred_network_dimensions=dimensions)
            flat.update(status='ok', **{key: summary[key] for key in
                        ['formula', 'density_g_cm3', 'cell_volume_angstrom3', 'atom_count']},
                        inferred_network_dimensions=';'.join(map(str, dimensions)))
            if include_porosity:
                geometry = pore_analysis(entry['structure'], probe_radius=probe_radius,
                                         samples=samples)
                row['porosity'] = geometry
                flat.update({key: geometry[key] for key in fields if key in geometry})
        except (ValueError, TypeError, KeyError, RuntimeError) as exc:
            row.update(status='error', error=str(exc))
            flat.update(status='error', error=str(exc))
        results.append(row)
        writer.writerow(flat)
    return {'kind': 'periodic_structure_comparison', 'entries': results,
            'successful': sum(row['status'] == 'ok' for row in results),
            'parameters': {'include_porosity': include_porosity,
                           'probe_radius_angstrom': probe_radius, 'samples': samples},
            'warnings': ['Comparisons use listed atoms and supplied cells. Model quality, '
                         'missing solvent/hydrogen, radii and sampling affect results.',
                         'Network dimensionality is distance-inferred; it does not assign an RCSR net.',
                         'Element site density counts atoms, not chemically accessible active sites.'],
            'exports': {'csv': table.getvalue()}}
